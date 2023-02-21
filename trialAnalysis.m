classdef trialAnalysis < handle 
    % Data intial processing for Joe Data Analysis
    % Filename:	trialAnalysis.m
    % Author:  James Hermus
    % Date:     1 Feb 2019
    % Description:	When initalized this class imports subject data and 
    % computes trial spesifict measures.
    
    properties
        
        % Subject Spesific Parameters
        subject % subject number in study
        trial % Trial Number (Feedback only) % ADD FIX THIS
        speed % Speed and direction are incorporated , according to speedNames
        names % subject idenifier
        speedNames % Name of speed + direction conditions
        l1a % Length of upper arm (m)
        l2a % Length of forearm (from elbow to center of fist) (m)
        ch  % Length from wrist to center of fist (m)
        wt  % Subject weight (kg)
        c1a % Distance from shoulder to uppre arm center of gravity (m)
        c2a % Distance from elbow to forearm center of gravity (m)
        I1a % Interia of upper arm (kg-m^2)
        I2a % Inertia of forearm  (including hand) (kg-m^2)
        m1a % Mass of upper arm (m)
        m2a % Mass of forearm (m)
        
        % Trial Spesific Parameters
        lc    % Crank Radius (m)
        I_crank % Inertia of the crank about pivot (kg-m^2)
        d1    % Horiz distance from crank center to shoulder center (m) (right when facing robot is positive)
        d2    % Vertical distance from crank center to shoulder center (m) (postive is toward the robot when facing it)
        B_crank
        
        % Assumed Human dynamics
        K % Stiffness matrix (N-m/rad)
        B % Damping matrix (N-m-s/rad)
        Z_gain % [G and lambda]  [2x1]
        turnSpeed
        
        sfrq % Sampling frequency
        
        % Measures
        t % Time (s)
        X % handel position: [Nx2] (m) in cartesian crank coordinates
        thcp  % Crank angular postion (rad)
        thcv  % Crank angular velocity (rad/s)
        F_rot % Force normal and tangental [Nx2] (m)
        q     % shoulder and elbow realitive joint angles [Nx2] (rad)
        q_dot % shoulder and elbow realitive joint velocities [Nx2] (rad)
        
        
        % Plotting params (add comments)
        armPosVec_X
        armPosVec_Y
        crank_x
        crank_y
        constraint_x
        constraint_y
        
        % Dependent Measures
        q_0   % ZFT postion in joint space [Nx2] (rad)
        q_0_dot
        X_0   % carteasion ZFT (m)
        Vt_0
        thcp_0
        X_0_mean
        lnrCoordinates
                
    end
    
    methods
        function  [this] = trialAnalysis(subject, speed, trial, Z_gain_input)
            this.subject = subject;
            this.speed = speed;
            this.trial = trial;
            if(nargin <= 3)
                % get Default Z in get_impedance4ZFTComp()
            elseif(nargin==4)
                this.Z_gain = Z_gain_input;
            elseif(nargin >=5)
                error('Too many inputs.');
            end
            
            % During construction always run
            this.sfrq = 200; % Enforeced sampling frequency
            this.get_subjectParam(); % Defines subject paramters
            this.get_subjectData();
            this.get_ZFT();
            this.get_curvatureVelData();
            this.get_lnrCoordinates();
            
        end
        
        %% Main analysis functions:
        function get_subjectParam(this)
            
            % General Shared Variables
            this.names = {'asn','bgn','iwd','onz','iuz','ktn','opo','inz','qba','idg'};
            this.speedNames = {'CW, Slow', 'CW, Med', 'CW, Fast', 'CCW, Slow', 'CCW, Med', 'CCW, Fast'};
            
            turnSpeedArray = [0.075,0.5,2,0.075,0.5,2];
            this.turnSpeed = turnSpeedArray(this.speed);
            
            %% HUMAN:
            [param, grpordr] = this.getparam();
            
            this.l1a = param(1);         % Length of upper arm
            this.l2a = param(2);         % Length of forearm (from elbow to center of fist)
            this.ch = param(3);          % From wrist to center of fist
            this.wt = param(4);          % Weight in kg
            this.lc = param(5);          % length of crank
            this.d1 = param(6);          % Horiz distance from shoulder to crank center
            this.d2 = param(7);          % Vert distance from shoulder to crank center
            
            %% CRANK:
            % For this application we desire the crank inertia with respect to the pivot not the center of mass
            this.I_crank = 3.716e-3; % Inertia of crank, measured (outboard subtracted).
            
            ka = 0.322*this.l1a;
            mf2 = 0.016*this.wt;
            cf2 = 0.430*(this.l2a-this.ch); kf2 = 0.303*(this.l2a-this.ch);
            mh = 0.006*this.wt + 0.5063; % Hand + (handle + outboard FT mass)
            
            this.c1a = 0.436*this.l1a; % Distance from shoulder to upper arm cg
            this.c2a = (mf2*cf2 + mh*this.l2a)/(mf2+mh); % Distance from elbow to forearm cg
            
            this.m1a = 0.028*this.wt; % mass of upper arm
            this.I1a = this.m1a*ka^2; % Inertia of upper arm
            
            If2 = mf2*kf2^2; % Inertia of forearm (minus hand)
            
            % Account for distance between hand mass and center of mass not being full forarm length.
            this.I2a = If2 + mh*(this.c2a-this.l2a)^2 + mf2*(this.c2a - cf2)^2; % Inertia of forearm (including hand)
            
            this.m2a = mf2+mh; % mass of forearm
            
            this.B_crank = 0;
            
            % Define joint torques using equlibrium trajectory
            this.get_impedance4ZFTComp();
            
            
        end
        
        function [] = get_impedance4ZFTComp(this)
            
            K11 = 29.5;
            K12 = 14.3;
            K22 = 39.3;
            
            if (~isempty(this.Z_gain))
                % Dont change K and B defined above
            else
                speed_StiffScale = [0.5,0.5,2,0.5,0.5,2]; % Scale for stiff conditions
                speed_DampingScale = [0.05,0.05,0.1,0.05,0.05,0.1]; % Scale for stiff conditions
                
                this.Z_gain = [speed_StiffScale(this.speed), speed_DampingScale(this.speed)];
            end
            
            G = this.Z_gain(1);
            lambda = this.Z_gain(2);
            
            this.K = G*[K11, K12; K12, K22];
            this.B = lambda*this.K;
        end
        
        function [] = get_subjectData(this)
            
            [param, grpordr] = this.getparam();
            jdex = find(this.speed == grpordr);
            
            % Select trial number skipping catch trials
            % Form Joe's code the blind matrix defines which trials did not
            % have visual feedback and which did. Assumes same pattern of
            % feedback for same speed and direction.
            blind = [0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0;...
                0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0;...
                0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0;...
                0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0;...
                0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0;...
                0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0];
            
            dex_Trial_Vision = find(~blind(this.speed,:));
            
            % Load data using:
            % jdex to switch to direction speed/direction
            % dex_Trial_Vision to choose only trials with had visual feedback
            % one is added to this.trial to skip the first trial in all cases.
            data = load(['crank_data/',...
                this.names{this.subject},'/', this.names{this.subject}, int2str(jdex), sprintf('%02d',dex_Trial_Vision(this.trial + 1)) ,'.ASC']);
            
            % Compute position as positive angle
            p = atan2(sin(data(:,1)), cos(data(:,1)));
            indxs = find(p(:,1)<0);
            p(indxs,1) = p(indxs) + 2*pi;
            
            % Define Fitler Motion
            cf = 10; % Cutoff frequency
            [b,a] = butter(2,cf/(this.sfrq/2));
            
            % Define Data position, velocity, and acceleration
            thcp = p; % Crank Angular Position
            thcv = filtfilt(b,a, data(:,2)*2*pi ); % Crank Angular Velocity 
                        
            % Define Fitler Normal Forces
            [b,a] = butter(2,cf/(this.sfrq/2));
            F_raw_Norm = data(:,3);
            F_raw_Norm_filt = filtfilt(b,a,F_raw_Norm);
            
            % Define Fitler Tangential Forces
            Tang_cf = [0.5,10,10,0.5,10,10];
            cf = Tang_cf(this.speed); % Cutoff frequency
            [b,a] = butter(2,cf/(this.sfrq/2));
            F_raw_Tang = data(:,4);
            F_raw_Tang_filt = filtfilt(b,a,F_raw_Tang);
            
            [ xe, ye, q1p, q2p ] = inverseKinomatics(this,thcp);
            t = [0:1/this.sfrq:1/this.sfrq*length(thcp)-1/this.sfrq]';
            
            %% Cut all data to only full cycles
            dexReset = find(abs(diff(thcp)) > 4);
            dexRange = [dexReset(1):dexReset(end)];
            
            this.thcp = thcp(dexRange);
            this.thcv = thcv(dexRange);
            this.F_rot = [F_raw_Norm_filt(dexRange)';F_raw_Tang_filt(dexRange)']; % Normal and Tangential Forces
            this.t = t(dexRange) - t(dexRange(1));
            this.X = [xe(dexRange), ye(dexRange)];
            this.q = [q1p(dexRange), q2p(dexRange)];
            
            %             figure;
            %             subplot(3,1,1); plot(this.t, this.thcp); ylabel('Postion (rad)');
            %             subplot(3,1,2); plot(this.t, this.thcv); ylabel('Velocity (rev/s)');
            %             subplot(3,1,3); plot(this.t, this.F_rot); ylabel('Force (N)');
            %             xlabel('Time (sec)');
            
            
        end
        
        %% Compute ZFT analysis
        function [] = get_ZFT(this)
            
            % Filename:	get_ZFT.m
            % Author:  James Hermus
            % Date:     9 June 2017
            % Description:	Computes the zero-force trajectory from subject
            % data with realitive damping.
            
            % Populate Plot feilds
            this.crank_x = [this.d1*ones(size(this.X(:,1))), this.d1 + this.lc*cos(this.thcp)];
            this.crank_y = [this.d2*ones(size(this.X(:,1))), this.d2 + this.lc*sin(this.thcp)];
            this.constraint_x = this.d1 - this.lc*cos(0:0.01:2*pi)';
            this.constraint_y = this.d2 - this.lc*sin(0:0.01:2*pi)';
            [this.armPosVec_X, this.armPosVec_Y] = plotPosArm(this,this.thcp);
            
            % Amply inverse kinomatics to get initial conditions
            [ xtmp, tmp, q1p_Start, q2p_Start, tmp ] = inverseKinomatics(this, 0 );
            
            x_init = [q1p_Start; q2p_Start];
            
            % Use ode45 to solve state equations
            [tvec, this.q_0] = ode45(@(t,y)this.getStatesZFT(t,y), this.t, x_init);
            
            % ZFT Postion
            [x_0, y_0] = forwKino(this,this.q_0(:,1),this.q_0(:,2));
            this.X_0 = [x_0, y_0];
            
            % ZFT Velocity
            for ii = 1:length(this.t)
                this.q_0_dot(ii,:) = getStatesZFT(this,this.t(ii),this.q_0(ii,:));
            end
            
        end
        
        function [] = get_curvatureVelData(this)
            % Filename:	get_curvatureVelData.m
            % Author:  James Hermus
            % Date:		May 18 2017
            % Description:	Compute vecloty consitent with Hermus et al. 
            % 2022 speed curvature paper.
            
            % Note for this entire function everything is in terms of X_0
            
            %% Compute and check for correctioness of higher dirivatives
            % Use Dohrmann spline diffrenetiation appraoch
            
            lo_bound = -20;
            hi_bound = 10;
            N_test = 35;
            
            B_vec = [10^-3, 10^-7, 10^-11, 10^-3, 10^-7, 10^-11]; 
            
            dorhman_B = B_vec(this.speed);
            [x,dx,ddx,dddx,Vx] = dohrmann(this.X_0(:,1), this.sfrq, dorhman_B);
            [y,dy,ddy,dddy,Vy] = dohrmann(this.X_0(:,2), this.sfrq, dorhman_B);
                        
            % Compter tangential velcoity
            vel = sqrt((dx).^2+(dy).^2);
            this.Vt_0 = vel;

            thcp_0 = atan2(this.X_0(:,2)-this.d2, this.X_0(:,1)-this.d1);
            indxs = find(thcp_0(:,1)<0);
            thcp_0(indxs,1) = thcp_0(indxs) + 2*pi;
            this.thcp_0 = thcp_0;
                        
        end
                
        function [] = get_lnrCoordinates(this)
                        
            this.X_0_mean = [mean(this.X_0(:,1)), mean(this.X_0(:,2))];

            % Subtract mean
            x = this.X_0(:,1) - this.X_0_mean(1);
            y = this.X_0(:,2) - this.X_0_mean(2);
                        
            % Compute lx,ly,lx',ly'
            for i = 1:4
                theta_d = (i-1)*(pi/4);
                l(i) = this.get_l(theta_d,x,y);
            end
            
            lx = l(1);
            ly = l(3);
            
            lx_prime = l(2);
            ly_prime = l(4);
            
            ln0 = log(lx/ly);
            ln45 = log(lx_prime/ly_prime);
            
            this.lnrCoordinates = [ln0;ln45];
            
        end
        
        function [l] = get_l(this,theta_d,x,y)
            
            % Assume mean zero
            if((size(x,2) ~= 1) || (size(y,2) ~= 1))
                error('transpose problem with get_l input');
            end
            
            % Computer theta
            theta_hat = atan2(y, x);
            indxs = find(theta_hat(:,1)<0);
            theta_hat(indxs,1) = theta_hat(indxs) + 2*pi;
            
            delta_theta = (2*pi)/200;
            dex0 = find( abs(wrapToPi(theta_hat-theta_d     )) < delta_theta );
            dex1 = find( abs(wrapToPi(theta_hat-(theta_d+pi))) < delta_theta );
            
            if(isempty(dex0) || isempty(dex1))
               % empty bin
               [tmp,dex0] = min(abs(theta_hat-theta_d));
               [tmp,dex1] = min(abs(theta_hat-(theta_d+pi)));

            end
            
            l = sqrt( ( mean(x(dex1)) - mean(x(dex0)) )^2 + ( mean(y(dex1)) - mean(y(dex0)) )^2 );
            
            % Check angle for atan2 and theta_d
            %             figure;plot(theta_hat);
            %
%                         figure;
%                         plot(x(dex0),y(dex0),'.r','markersize',25); hold on;
%                         plot(x(dex1),y(dex1),'.b','markersize',25); hold on;
%                         plot(x,y,'.'); hold on;
            
        end
                
        %% Helper Functions:
        function qdv = getStatesZFT(this,tt,xx)
            % Forward dynamic model of a two-link planar mechanism
            % Neville Hogan, 2007
            % Based on a script by Jerome J. Palazzolo, 2003
            
            % % Extract states
            q1dp=xx(1);
            q2dp=xx(2);
            
            % Follow Defined q1p, q1v, q1a / q2p, q2v, q2a / thcp, thcv, thca
            thcp = interp1(this.t, this.thcp, tt);
            thcv =  interp1(this.t, this.thcv, tt);
            F_rot(1,:) = interp1(this.t, this.F_rot(1,:), tt);
            F_rot(2,:) = interp1(this.t, this.F_rot(2,:), tt);
            
            % Amply inverse kinomatics
            [ xe, ye, q1p, q2p, thea ] = inverseKinomatics(this, thcp);
            
            % Define Jacobian of the arm
            [J] = get_jacobian(this,[q1p,q2p]);
            
            % nn % Normal unit vector
            % ee % tangential unit vector
            ee = [ - sin(thcp); cos(thcp) ];
            nn = [  cos(thcp); sin(thcp) ];
            
            qv = inv(J)*(this.lc*thcv*ee);
            q1v = qv(1);
            q2v = qv(2);
            
            % Define inertial matrix of the arm
            [M] = get_massMatrixArm(this, q1p, q2p);
            
            % Define centifugal and Coriolis forces of the arm
            h = [ -this.m2a*this.l1a*this.c2a*sin(q2p) * (2*q1v*q2v + q2v^2);...
                this.m2a*this.l1a*this.c2a*sin(q2p) * q1v^2];
            
            % Define Jacobian Dot of the arm
            Jd11 = -( this.l1a*q1v*cos(q1p) + this.l2a*(q1v+q2v)*cos(q1p + q2p) );
            Jd12 = -this.l2a*(q1v+q2v)*cos(q1p + q2p);
            Jd21 = -this.l1a*q1v*sin(q1p) - this.l2a*(q1v+q2v)*sin(q1p + q2p);
            Jd22 = -this.l2a*(q1v+q2v)*sin(q1p + q2p);
            Jdot = [Jd11, Jd12; Jd21, Jd22];
            
            % configuration-dependent damping term
            H = this.B_crank*this.thcv + this.lc * ee'* inv(J') * (h - M * inv(J) * ( this.lc * thcv^2 * nn + Jdot * [q1v; q2v] ) );
            
            % Force implimentation
            F_cart = [F_rot(1)*nn + F_rot(2)*ee];
            C1 = M*inv(J) * ( ( J*inv(M)*J' + this.lc^2*inv(this.I_crank)*ee*ee' )*F_cart - Jdot*[q1v;q2v] - this.lc*thcv*(thcv*nn + this.B_crank*inv(this.I_crank)*ee) ) + h;
            qdv = inv(this.B)*(C1 - this.K*([q1dp;q2dp] - [q1p;q2p]) ) + [q1v;q2v];
        end
                
        function [tau_sys, deltaQ, deltaQ_dot] = getTorque4K(this,ii)
            
            q1dp = this.q_0(ii,1);
            q2dp = this.q_0(ii,2);
            q1dv = this.q_0_dot(ii,1);
            q2dv = this.q_0_dot(ii,2);
            
            thcp = this.thcp(ii);
            thcv =  this.thcv(ii);
            F_rot(1,:) = this.F_rot(1,ii);
            F_rot(2,:) = this.F_rot(2,ii);
            
            % Amply inverse kinomatics
            [ xe, ye, q1p, q2p, thea ] = inverseKinomatics(this, thcp);
            
            % Define Jacobian of the arm
            [J] = get_jacobian(this,[q1p,q2p]);
            
            % nn % Normal unit vector
            % ee % tangential unit vector
            ee = [ - sin(thcp); cos(thcp) ];
            nn = [  cos(thcp); sin(thcp) ];
            
            qv = inv(J)*(this.lc*thcv*ee);
            q1v = qv(1);
            q2v = qv(2);
            
            % Define inertial matrix of the arm
            [M] = get_massMatrixArm(this, q1p, q2p);
            
            % Define centifugal and Coriolis forces of the arm
            h = [ -this.m2a*this.l1a*this.c2a*sin(q2p) * (2*q1v*q2v + q2v^2);...
                this.m2a*this.l1a*this.c2a*sin(q2p) * q1v^2];
            
            % Define Jacobian Dot of the arm
            Jd11 = -( this.l1a*q1v*cos(q1p) + this.l2a*(q1v+q2v)*cos(q1p + q2p) );
            Jd12 = -this.l2a*(q1v+q2v)*cos(q1p + q2p);
            Jd21 = -this.l1a*q1v*sin(q1p) - this.l2a*(q1v+q2v)*sin(q1p + q2p);
            Jd22 = -this.l2a*(q1v+q2v)*sin(q1p + q2p);
            Jdot = [Jd11, Jd12; Jd21, Jd22];
            
            % configuration-dependent damping term
            H = this.B_crank*this.thcv + this.lc * ee'* inv(J') * (h - M * inv(J) * ( this.lc * thcv^2 * nn + Jdot * [q1v; q2v] ) );
            
            % Force implimentation
            F_cart = [F_rot(1)*nn + F_rot(2)*ee];
            
            tau_sys = M*inv(J) * ( ( J*inv(M)*J' + this.lc^2*inv(this.I_crank)*ee*ee' )*F_cart - Jdot*[q1v;q2v] - this.lc*thcv*(thcv*nn + this.B_crank*inv(this.I_crank)*ee) ) + h;
            
            deltaQ = ([q1dp;q2dp] - [q1p;q2p]);
            deltaQ_dot = ([q1dv;q2dv] - [q1v;q2v]);
            
        end
        
        function [x, y] = forwKino(this,q1p,q2p)
            
            % Requires postion data in shoulder coordinates
            % Exports forward kinomatics in Shoulder centered coordinates
            
            x = [this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p)];
            y = [this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p)];
            
        end
        
        function [ xe, ye, q1p, q2p, thea ] = inverseKinomatics( this, thc )
            %Compute joint angles
            %   Computes crank postion and joint angles input params provides the
            %   subject data from the get"kino number"param()
            
            % Compute crank postion in carteasion
            xe = this.d1 + this.lc*cos(thc);
            ye = this.d2 + this.lc*sin(thc);
            ld = sqrt(xe.^2 + ye.^2);
            
            % Right Hand
            alpha1 = atan2(ye,xe);
            alpha2 = acos((xe.^2 + ye.^2 + this.l1a^2 - this.l2a^2)./(2*this.l1a*sqrt(xe.^2 + ye.^2)));
            alpha3 = acos((xe.^2 + ye.^2 + this.l2a^2 - this.l1a^2)./(2*this.l2a*sqrt(xe.^2 + ye.^2)));
            
            % spatial angles
            q1p = alpha1 - alpha2;
            q2p = alpha1 + alpha3 - q1p;
            
            %% Old Code that works
            % % Compute elbow angle relative
            % q2p = pi - acos((-ld.^2+l1a^2+l2a^2)./(2*l1a*l2a));
            %
            % % shoulder angle
            % q1p =  atan2(ye,xe) - acos((ld.^2+l1a^2-l2a^2)./(2.*ld.*l1a));
            
            % elbow angle absolute
            thea = q1p + q2p;
            
        end
        
        function [J] = get_jacobian(this,q)
            
            q1p = q(1);
            q2p = q(2);
            
            J11 = -( this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p) );
            J12 = -this.l2a*sin(q1p + q2p);
            J21 = this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p);
            J22 = this.l2a*cos(q1p + q2p);
            J = [J11, J12; J21, J22];
            
        end
        
        function [M] = get_massMatrixArm(this, q1p, q2p)
            
            % Define inertial matrix of the arm
            M11 = this.m1a*this.c1a^2 + this.m2a*( this.l1a^2 + this.c2a^2 + 2*this.l1a*this.c2a*cos(q2p) ) + this.I1a + this.I2a;
            M12 = this.m2a*( this.c2a^2 + this.l1a*this.c2a*cos(q2p) ) + this.I2a;
            M22 = this.m2a*this.c2a^2 + this.I2a;
            M = [ M11, M12; M12, M22 ];
            
        end
        
        function [armPosVec_X, armPosVec_Y] = plotPosArm(this, thcp)
            
            % Takes in crank coorinates
            [ xe, ye, q1p, q2p, thea ] = inverseKinomatics( this, thcp );
            
            % Returns the arm postion for plotting
            armPosVec_X = [zeros(length(q1p),1), this.l1a*cos(q1p),...
                this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p)];
            
            armPosVec_Y = [zeros(length(q1p),1), this.l1a*sin(q1p),...
                this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p)];
            
            %             figure; plot(armPosVec_X, armPosVec_Y,'-ok'); axis equal; grid on;
            
        end
        
        function [V,D] = get_eigs_OrderByMagnitude(this,A)
            
            [V,D] = eig(A);
            D = diag(D);
            [D,dex] = sort(abs(D),'descend'); % Reset D order
            V = [V(:,dex(1)), V(:,dex(2))];   % Reset V order
            
        end
        
        function [] = get_zftPlot(this)
            % Zero Force Trajectory Plot
            figure('position',[693 279 560 420]);
            pointsize = 50;
            LW = 3.5;
            FS = 18;
            dexRange = 100:length(this.X_0);
            %             if(this.speed == 1 || this.speed == 4)
            %                 plot(this.X_0_static(:,1)-this.d1,this.X_0_static(:,2)-this.d2,'-r','linewidth',2.5); hold on;
            %             end
            scatter(this.X_0(dexRange,1)-this.d1, this.X_0(dexRange,2)-this.d2, pointsize,(this.Vt_0(dexRange)-min(this.Vt_0(dexRange)))/(max(this.Vt_0(dexRange))-min(this.Vt_0(dexRange))),'filled'); hold on;
            hcb = colorbar;
            plot(this.lc*cos(0:0.05:2*pi),this.lc*sin(0:0.05:2*pi),'--k','linewidth',LW);hold on;
            set(gca,'FontSize',FS);
            
            ylabel(hcb,'V(s) (m/s)','fontsize',FS); ylabel('Y-Position (m)'); xlabel('X-Position (m)');
            ylabel(hcb,'Normalized Speed');
            ylabel('Y-Position (m)');
            xlabel('X-Position (m)');axis equal; %grid on;
            ylim([-0.18 0.18]);
            xlim([-0.18 0.18]);
            set(gca,'FontSize',28);
            caxis([0 1]);
            
            %             figure;plot(this.Vt_0(dexRange));
            
            %             figure('position',[114 279 560 420]);
            %             subplot(2,1,1); plot(this.thcp,(this.lc - sqrt((this.X_0(:,1)-this.d1).^2 + (this.X_0(:,2)-this.d2).^2))*(1000),'.');
            %             ylabel('\Delta x radial (mm)'); grid on;xlim([0 2*pi]); set(gca,'fontsize',18);
            %
            %             subplot(2,1,2); plot(this.thcp,this.F_rot(1,:),'.');
            %             xlabel('\theta_c (rad)'); ylabel('Normal Force (N)'); grid on; xlim([0 2*pi]); set(gca,'fontsize',18);
            
            
        end
        
        function [] = save_zftAnimation(this,videoFileName)
            
            FS = 18;
            figure;
            v = VideoWriter(videoFileName, 'MPEG-4');
            skipSamples = 10;
            v.FrameRate = this.sfrq/skipSamples;
            open(v);
            
            for i = 1:skipSamples:length(this.t)
                
                ax = plot(this.armPosVec_X(i,:), this.armPosVec_Y(i,:), '-b',...                    
                    this.crank_x(i,:),     this.crank_y(i,:),     '-k',...
                    this.constraint_x,     this.constraint_y,     ':k',...
                    this.armPosVec_X(i,:), this.armPosVec_Y(i,:), '.k',...
                    this.crank_x(i,:),     this.crank_y(i,:),     '.k',...
                    this.X_0(i,1),           this.X_0(i,2),           '.r','linewidth',5,'MarkerSize',50);
                axis equal;
                ylim([0 0.8]); xlim([-0.4 0.4]);
                ylabel('Y Distance (m)','fontsize',14);
                xlabel('X Distance (m)','fontsize',14);
                ax(1).LineWidth = 7;
                ax(3).LineWidth = 2.5;
                ax(6).MarkerSize = 40;
                legend([ax(1),ax(2),ax(3),ax(6)],{'Arm','Crank','Constraint Path','Zero-Force Trajectory'},'location','northwest','fontsize',FS-4);
                    
                axis equal;%grid on; 
                ylim([0 0.8]); xlim([-0.4 0.4]);
                ylabel('Y Distance (m)','fontsize',FS);
                xlabel('X Distance (m)','fontsize',FS);
                set(gca,'box', 'off','linewidth',2,'Fontsize',FS);
                set(gcf,'Color',[1,1,1]);
%                 legend('Arm','Zero-Force Trajectory','Crank','Constraint Path','location','northwest');

%                 set(gca,'YTickLabel',[],'XTickLabel',[]);


                frame = getframe(gcf);
                writeVideo(v,frame);
            end
            close(v);
            
        end
       
        function [ X_bin, vel_bin ] = customBinning(this, x, y, vel, p)
            
            N = 200; % number of bins
            posEdges = linspace(0,2*pi,N+1);
            
            for j = 1:N
                
                binDex = find(p > posEdges(j) & p <= posEdges(j+1)); % Find index pos prosition with in each bin
                inda(j,1:length(binDex)) = binDex;                   % save index ast inda
                numPerBin(j) = length(binDex);
                vel_bin(j) = mean(vel(binDex));                  % take mean of velocity positon points with in each bin
                vel_bin_tot{j} = vel(binDex);
                
                % Find mean radius too
                x_bin(j) = mean(x(binDex));
                y_bin(j) = mean(y(binDex));
                r(j) = sqrt(mean(x(binDex)).^2+mean(y(binDex)).^2);
                
            end
            X_bin = [x_bin; y_bin];
            
        end
        
        %% Get paramerters for spesific subjects
        % Recycle structure from Joe
        function [param, grpordr] = getparam(this)
            
            % function [param, grpordr] = getparam(subj)
            %
            % Another way of gettin subject parameters -- like a front end to
            % the different get'subj'params.
            
            if(    this.subject == 1 )
                [param, grpordr] = this.getasnparam();
            elseif(this.subject == 2 )
                [param, grpordr] = this.getbgnparam();
            elseif(this.subject == 3 )
                [param, grpordr] = this.getiwdparam();
            elseif(this.subject == 4 )
                [param, grpordr] = this.getonzparam();
            elseif(this.subject == 5 )
                [param, grpordr] = this.getiuzparam();
            elseif(this.subject == 6 )
                [param, grpordr] = this.getktnparam();
            elseif(this.subject == 7 )
                [param, grpordr] = this.getopoparam();
            elseif(this.subject == 8 )
                [param, grpordr] = this.getinzparam();
            elseif(this.subject == 9 )
                [param, grpordr] = this.getqbaparam();
            elseif(this.subject == 10 )
                [param, grpordr] = this.getidgparam();
            else
                error('Subject get params not selected.');
            end
            
        end
        
        function [param, grpordr] = getasnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,3,4,6,5,1];
            
            la = 14.5*0.0254;   % Length of upper arm
            lf = 14.75*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 175*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 18.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getbgnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [6,4,2,5,1,3]; % Omit set '0'
            
            la = 14.25*0.0254;   % Length of upper arm
            lf = 14.5*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 18.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getidgparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,4,3,6,1,5]; % Omit set '0'
            
            la = 12.75*0.0254;   % Length of upper arm
            lf = 13.0*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 173*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.125*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getinzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,6,5,3,4,2]; % Omit set '0'
            
            la = 14.25*0.0254;   % Length of upper arm
            lf = 15.0*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 205*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.75*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getiuzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,4,3,6,5,1]; % Omit set '0'
            
            la = 13.75*0.0254;   % Length of upper arm
            lf = 14*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getiwdparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,3,5,4,6,2]; % Omit set '0'
            
            la = 11.5*0.0254;   % Length of upper arm
            lf = 13.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3*0.0254;   % From wrist to center of fist
            wt = 185*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 16.75*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getktnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [5,1,6,2,3,4]; % Omit set '0'
            
            la = 14.75*0.0254;   % Length of upper arm
            lf = 14.75*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 155*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
            
        end
        
        function [param, grpordr] = getonzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,3,1,4,5,6]; % Omit set '0'
            
            la = 13.5*0.0254;   % Length of upper arm
            lf = 14.5*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 155*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getopoparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,1,6,3,4,5]; % Omit set '0'
            
            la = 14.5*0.0254;   % Length of upper arm
            lf = 15.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
            
        end
        
        function [param, grpordr] = getqbaparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,6,5,3,2,4]; % Omit set '0'
            
            la = 13.5*0.0254;   % Length of upper arm
            lf = 13.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 180*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
    end
end