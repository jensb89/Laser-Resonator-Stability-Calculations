classdef Resonator < handle
    % Calculate Stability of a laser resonator with ABCD matrix formalism
    % 
    % Syntax:
    % Resonator(dK,nK,thetaK,R,L,lambda,pol)
    %
    % dK = (TiSa-) crystal thickness / m
    % nK = crystal refractive index
    % thetaK = angle of incidence (e.g. Brewster angle)
    % R = Radius of curvature (ROC) of the Cavity mirrors
    % L = Length of the resonator
    % lambda = Wavelength
    % pol = Polarisation ('s' for sagittal plane or 'p' for tangential plane)
    %
    % Version 1.0
    %
    % Copyright Jens Brauer, 03.2014
    
    properties
        nK = 1.76
        dK = 3E-3
        thetaK =   1.0541     %angle of incidence on crystal (Brewster) / rad
        R  =   0.1            %focusing mirrors radius of curvature R=2*f / m 
        L_ges  =   1.84       %resonator length / m (L_arm1+ L_arm2 + R + dK
        
        %Length of Arms
        L_arms %Larm1+Larm2
        L1 %Length of Arm1
        L2 %Length of Arm2
        
        %Folding angles
        theta_ges
        theta1
        theta2
        
        % Cavity dimensions
        Delta = 20E-4 %deviation of distance between cavity mirrors from 
                      %confocal configuration (2*f+delta) / m
        dSK1 %Length from Cavity Mirror 1 to Crystal
        dSK2 %Length from Cavity Mirror 2 to Crystal
        d_eff %Effective path in tilted crystal  %unused (d/cos(theta_in)) theta_in = asin(sin(thetaK)/nK)
        
        v_rep  %repitition rate / Hz
        
        lambda =   800e-9       %wavelength / m
        c  =   2.998e8          %speed of light in vaccum
        
        q1 % Eigenmode of the system (Imaginary number 1/q = 1/R - i lambda/(pi * w0^2))
        q2
        
        pol = 's' %Polarisation ('s' or 'p')
    end
    
    methods
        
        function this = Resonator(dK,nK,thetaK,R,L,lambda,pol)
            this.dK = dK;
            this.nK = nK;
            this.thetaK = thetaK;
            this.R = R;
            this.L_ges = L;
            this.lambda = lambda;
            this.pol = pol;
            
            this.L_arms = L-R-dK;
            this.v_rep = this.c/(2*this.L_ges);
            this.dSK1 = 0.5*(R - dK + this.Delta);
            this.dSK2 = this.dSK1;
        end
        
        function theta_ges_ = getTheta(this)
            %controlling the sum of folding angles in both resonator arms
            %find optimal angle for astigmatism free beam
            %it is : 
            %(n^2-1)*sqrt(n^2+1)/n^4 * dK = R sin(theta_ges/4) * %tan(theta_ges/4)
            %
            %Optimal angle solves the equation:
            asti = @(x) abs((this.nK^2-1)*sqrt(this.nK^2+1)/this.nK^4*this.dK/this.R...
                            -sin(pi/180*x/4)*tan(pi/180*x/4));
            theta_ges_    = fminbnd(asti,20,40)*pi/180;  
            this.theta_ges = theta_ges_;
        end
        
        function setArmLengthDist(this,val)
            % Set Arm Length Distribution
            if val>1
                val = val/100;
            end
            this.L1 = this.L_arms * val;
            this.L2 = this.L_arms * (1-val);
        end
        
        function setThetaDist(this,val)
            % Set Folding Angle Distribution
            if val>1
                val = val/100;
            end
            if isempty(this.theta_ges)
                this.getTheta();
            end
            
            this.theta1 = this.theta_ges * val; 
            this.theta2 = this.theta_ges * (1-val);
        end
        
        function setPolarisation(this, pol)
            this.pol = pol;
        end
        
        function setRelCrystalPos(this, relpos)
            this.dSK1 = relpos*(this.R-this.dK+this.Delta);
            this.dSK2 = (1-relpos)*(this.R-this.dK+this.Delta);
        end
            
        function setDelta(this, delta)
            %Set deviation of cavity mirror distance from 2f
            eta_actual = this.dSK1/(this.R-this.dK+this.Delta);
            
            this.Delta = delta;
            % Update distances to crystal as well:
            this.dSK1 = eta_actual*(this.R + this.Delta - this.dK);
            this.dSK2 = (1-eta_actual)*(this.R + this.Delta - this.dK);
        end
        
        function [M_arm1,M_arm2,w_arm1,w_arm2,S_arm1,S_arm2] = calcRoundTrip(this)
            % M1---L1----MC1----LdSK1---K---LdSK2---MC2----L2---M2
            % M = Mirror
            % MC = Curved Mirror from the Cavity
            % K = Crystal
            % L = Length of Arm
            
            M1 = mirror();
            M2 = mirror();
            MC1 = curved_mirror(this.R,this.theta1/2,this.pol); %/2 -> Halber Einfallswinkel (denn Strahl auf optischer Achse und nicht Theta_in = Theta_out wie per Defintiont)(vgl Kasper oder "Siegman, Lasers")
            MC2 = curved_mirror(this.R,this.theta2/2,this.pol);
            L1 = free(this.L1);
            L2 = free(this.L2);
            LdSK1 = free(this.dSK1);
            LdSK2 = free(this.dSK2);
            K = tilted_crystal(this.dK,this.thetaK,this.nK,this.pol);
            
            M_round1 = L1*MC1*LdSK1*K*LdSK2*MC2*L2;
            M_round2 = L2*MC2*LdSK2*K*LdSK1*MC1*L1;
            
            M_arm1 = M_round1 * M2 * M_round2;
            M_arm2 = M_round2 * M1 * M_round1;
            
            
            % Calculate Eigenmode of the resonator 
            % Arm1
            [A,B,C,D] = this.getABCD_Components(M_arm1);
            S_arm1 = abs((A+D)/2); %Stability Factor (Stable for S<1)
            [q_arm1,w_arm1,R_arm1] = calcResonatorEigenmode(this,A,B,D);
            
            % Arm2
            [A,B,C,D] = this.getABCD_Components(M_arm2);
            S_arm2 = abs((A+D)/2);
            [q_arm2,w_arm2,R_arm2] = calcResonatorEigenmode(this,A,B,D);
            
            this.q1 = q_arm1;
            this.q2 = q_arm2;
            
        end
        
        
        function [q,w,R] = calcResonatorEigenmode(this,A,B,D)
            % Returns 1/q for the Eigenmode   
            q = (D-A)/(2*B) + 1/B * sqrt(((A+D)/2)^2-1); %1/q
            q = 1/q;
            q = real(q) + 1i*abs(imag(q)); %abs because of quadratic eigenmode equation, solution +-
            R = real(q);
            w = sqrt(-this.lambda/pi/imag(1/q));
        end
        
        function estimation_Delta(this)
            %Calculate estimated stability gap
            f = this.R/2;
            delta1=f^2/(this.L1-f);
            delta2=f^2/(this.L2-f);
            delta_max = delta1+delta2;
            fprintf('Estimated stability : \n Stable from %1.2fmm to %1.2fmm and from %1.2fmm to %1.2fmm\n',...
                0,delta1/1E-3,delta2/1E-3,delta_max/1E-3);
        end
        
        
        function [z,w] = beamInsideCrystal(this)
            % ABCD Matrix Formalism from M2 to K (right before the crystal)
            % Then the Gaussian beam inside the crystal is calculated with
            % q2 = Aq1+B / (Cq1+D)
            L2 = free(this.L2);
            LdSK2 = free(this.dSK2);
            MC2 = curved_mirror(this.R,this.theta2/2,this.pol);
            if strcmp(this.pol,'p')
                nVak = 1;
                theta_out = this.thetaK;
                theta_in = asin(sin(theta_out)*nVak/this.nK);
                %K_in = eye(2);
                K_in = [cos(theta_out)/cos(theta_in) 0;0 cos(theta_in)/cos(theta_out)]; %Curved_Interface, arbitrary incidence, with R=inf
            else
                K_in = eye(2);
            end
            M = K_in * LdSK2 * MC2 * L2;%L2 * MC2 * LdSK2 * K_in;
            
            [A,B,C,D] = this.getABCD_Components(M);
            q0 = this.q2; %Eigenmode at M2
            z = linspace(0,2*this.dK,300);
            q1 = (A*q0+B)/(C*q0+D);
            
            w  = sqrt(-this.lambda/pi./imag(1./(q1+z/this.nK))); 
        end
        
        function [z_pump, w_pump] = PumpBeamInsideCrystal(this,lambda_p,w0,f1,f2,dCM,L12,L3)
            %---()----()-----|(----/ /
            %   f1 L12 f2 L3  f3 L4 K_in
            %
            % w0 = Beam diameter pump beam
            % f1 = Focus length lens1 (inf if not used)
            % f2 = Focus length lens2
            % dCM = Thickness of (curved) Cavity mirror
            % L12,L3,L4 = path length between elements
            % 
            % Output = z and w(z) beam diameter inside the crystal
            
            L1 = free(10E-2); %Before first lens
            L12 = free(L12);
            L3 = free(L3);
            L4 = free(this.dSK1);
            f1 = lens(f1);
            f2 = lens(f2);
            theta = this.theta2/180*pi;
            n = 1.5; % Cavity mirror index of refraction
            f3 = curved_mirror_transmission(this.R,dCM,theta,n,this.pol);

            if strcmp(this.pol,'p')
                nVak = 1;
                theta_out = this.thetaK;
                theta_in = asin(sin(theta_out)*nVak/this.nK);
                K_in = [cos(theta_out)/cos(theta_in) 0;0 cos(theta_in)/cos(theta_out)]; %Curved_Interface, arbitrary incidence, with R=inf
            else
                K_in = eye(2);
            end

            M = K_in * L4 * f3 * L3 * f2 * L12 * f1 * L1;
            
            z_pump=linspace(0,2*this.dK,300);

            [A,B,C,D] = this.getABCD_Components(M);
            q1 =  1/(-1i*lambda_p/pi/w0^2); %Incoming beam
            q2 = (A*q1+B)/(C*q1+D);

            nK=this.nK;
            w_pump  = sqrt(-lambda_p/pi./imag(1./(q2+z_pump/nK))); 
        end
        
    end
        
    
    methods(Static)
        function [A,B,C,D] = getABCD_Components(ABCD_Matrix)
            A = ABCD_Matrix(1,1);
            B = ABCD_Matrix(1,2);
            C = ABCD_Matrix(2,1);
            D = ABCD_Matrix(2,2);
        end
        
        function qfactor=q_(w,R,varargin)
            % Returns the q-factor of a Gaussian beam
            % w      = 1/e Field radius 
            % R      = Radius of curvature of phasefront
            % lambda = wavelength
            % Compare Kogelnik, 1966

            if nargin>=3, lambda=varargin{1}; else lambda=800E-9; end

            qfactor = 1./R - 1i.*lambda/(pi.*w.^2); %1/q
            qfactor = 1/qfactor;                    %q     
        end
    end
    
    
end

