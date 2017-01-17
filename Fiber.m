classdef Fiber < matlab.System
    % Fibre Add summary here
    %
    
    properties
        % Public, tunable properties.
        
        % wavelength of source in meter
        lambda; 
        % core radius in meter
        r; 
        % core refrative index
        %Nc = 1.456  1.5699;
        Nc;
        % cladding refractive index
        Ng;
        
        % Wave number
        k;
        % Cutoff Frequency
        V;
        % Index difference
        delta;
        % Eigenvalues of TM modes
        uTM; NeffTM;
        % Eigenvalues of TE modes
        uTE; NeffTE;
        % Eigenvalues of EH modes
        uEH; NeffEH;
        % Eigenvalues of HE modes
        uHE; NeffHE;
        % v and m max
        vmax; mmax;
        
    end
    
    properties (DiscreteState)
    end
    
    properties (Access = private)
        % Pre-computed constants.        
    end
    
    methods
       function obj = Fiber(lambda,r,Nc,Ng,vmax, mmax)
            if nargin < 5
                vmax = 3;mmax = 5;
                obj.vmax = vmax;
                obj.mmax = mmax;
            end
            obj.vmax = vmax;
            obj.mmax = mmax;

            format long
            obj.lambda = lambda;
            obj.r = r;
            obj.Nc = Nc;
            obj.Ng = Ng;
            obj.k = 2*pi/lambda;
            obj.V = obj.k*r*sqrt(Nc^2 - Ng^2);
            obj.delta = (Nc^2 - Ng^2)/(2*Nc^2);

            % Eigenvalues of EH modes
            obj.uEH = zeros(mmax,vmax); obj.NeffEH = zeros(mmax,vmax);
            % Eigenvalues of HE modes
            obj.uHE  = zeros(mmax,vmax); obj.NeffHE  = zeros(mmax,vmax);
        end
                
        function eigenValuesTE(obj)
            syms U;
            W = sym(@(U) sqrt(obj.V^2 - U^2));
            W1 = sym(@(U1) sqrt(obj.V^2 - U1^2));
            
            eq = sym(@(U1)( (besselj(1,U1)/(U1*besselj(0,U1))) + (besselk(1,W1)/(W1*besselk(0,W1))) ));
            f(U) = ( (besselj(1,U)/(U*besselj(0,U))) + (besselk(1,W)/(W*besselk(0,W))) );
            
            v = 0;
            sol = findZeros(obj,f,eq,v);
            [obj.uTE, obj.NeffTE] = setData(obj, sol);
        end
        
        function eigenValuesTM(obj)
            U = sym('U');
            W = sym(@(U) sqrt(obj.V^2 - U^2));
            W1 = sym(@(U1) sqrt(obj.V^2 - U1^2));
            
            eq = sym(@(U1) ((obj.Nc^2*besselj(1,U1))/(U1*besselj(0,U1))) + ((obj.Ng^2*besselk(1,W1))/(W1*besselk(0,W1))) );
            f(U) =  ((obj.Nc^2*besselj(1,U))/(U*besselj(0,U))) + ((obj.Ng^2*besselk(1,W))/(W*besselk(0,W))) ;
            
            v = 0;     
            sol = findZeros(obj,f,eq,v);
            [obj.uTM, obj.NeffTM] = setData(obj, sol);
        end
        
        function eigenValuesHyb(obj)          
            syms U v;
            W = sym(@(U) sqrt(obj.V^2 - U^2));
            W1 = sym(@(U1) sqrt(obj.V^2 - U1^2));
            
            beta = sym(@(U) sqrt(obj.k^2*obj.Nc^2 - U^2/obj.r^2));
            beta1 = sym(@(U1) sqrt(obj.k^2*obj.Nc^2 - U1^2/obj.r^2));
            
            neff = sym(@(U) beta/obj.k);
            neff1 = sym(@(U1) beta1/obj.k);
            
            x = sym(@(U, v) (diff(besselj(v,U),U))/(U*besselj(v,U)));
            x1 = sym(@(U1, v) (diff(besselj(v,U1),U1))/(U1*besselj(v,U1)));
            
            c = sym(@(U, v) ((v*neff)^2)*((obj.V/(U*W))^4));
            c1 = sym(@(U1, v) ((v*neff1)^2)*((obj.V/(U1*W1))^4));
            
            b = sym(@(U, v) ((v*besselk(v, W))/W - besselk(v + 1, W))/(W*besselk(v,W)));
            b1 = sym(@(U1, v) ((v*besselk(v, W1))/W1 - besselk(v + 1, W1))/(W1*besselk(v,W1)));            
            
            fhe(U) = -b*(1-obj.delta) - sqrt(b^2*obj.delta^2 + c/obj.Nc^2) - x;
            feh(U) = -b*(1-obj.delta) + sqrt(b^2*obj.delta^2 + c/obj.Nc^2) - x;                         
            
            for v = 1:1:obj.vmax
                eqhe = sym(@(U1) eval(-b1*(1-obj.delta) - sqrt(b1^2*obj.delta^2 + c1/obj.Nc^2) - x1));
                eqeh = sym(@(U1) eval(-b1*(1-obj.delta) + sqrt(b1^2*obj.delta^2 + c1/obj.Nc^2) - x1));

                sol = findZeros(obj,fhe,eqhe,v);
                [obj.uHE(:,v), obj.NeffHE(:,v)] = setData(obj, sol);
                sol = findZeros(obj,feh,eqeh,v);
                [obj.uEH(:,v), obj.NeffEH(:,v)] = setData(obj, sol);
            end
        end

        
        function [u, Neff] = setData(obj, sol)
            if numel(sol) > obj.mmax
                error('Value of mmax must be at least %d', numel(sol))
            end
            beta = sym(@(U) sqrt(obj.k^2*obj.Nc^2 - U^2/obj.r^2));
            neff = sym(@(U) beta/obj.k);
            
            u = zeros(obj.mmax,1);
            Neff = zeros(obj.mmax,1);
            for soli = 1:1:numel(sol)
               u(soli,1) = sol(soli);
               U = sol(soli); %#ok<NASGU>
               Neff(soli,1) = eval(neff);
            end

        end
        
        function sol = findZeros(obj,f,eq,v)         
            sol = zeros(0,10);
            soli = 0;
            err = 1e-2;
            pas = 1e-1;
            U1 = sym('U1');
            tmp2 = pi;
            
            for U = pas:pas:obj.V
                if( abs(U*besselj(v,U)) >= err)
                    tmp1 = tmp2;
                    tmp2 = eval(f(U));      
                    if(tmp1 ~= pi && sign(tmp1) ~= sign(tmp2))  
                        soli = soli + 1;
                        sol(soli) = vpasolve(eq == 0, U1, U);          
                    end            
                end
            end            
            sol = unique(sol);
        end
    end
end
