function [A,B, z] = edogral(T,Tc,P,Pc,w,ee)

Tr = T/Tc;
Pr = P/Pc;

switch ee
    case 1
        A = (27/64)*Pr/(Tr^2)
        B = Pr/(8*Tr)
        p = -1-B
        q = A
        r = -A*B
        z = roots([1 p q r])
        z = z(imag(z)==0)
        
        uu = -A/z;
        hh = z-1-A/z;
        ss = ln(z-B);
        theta = A/z;
        phi = exp(z-1-log(z-B)-theta)
        
        disp('H-H* / RT = ')
        disp(hh)
        disp('U-U* / RT = ')
        disp(uu)
        disp('S-S* / R = ')
        disp(ss)
        disp('Θ= ')
        disp(theta)
        disp('φ= ')
        disp(phi)
        
    case 2
        A = 0.42748*(Pr/Tr^(2.5))
        B = 0.08664*(Pr/Tr)
        p = -1
        q = A-B-B^2
        r = -A*B
        z = roots([1 p q r])
        z = z(imag(z)==0)
        
        uu = -3*A/(2*B)*log(1+B/z);
        hh = z-1+uu;
        ss = log(z-B)-A/(2*B)*log(1+B/z);
        theta = (A/B)*log(1+B/z);
        phi = exp(z-1-log(z-B)-theta)
        
        disp('H-H* / RT = ')
        disp(hh)
        disp('U-U* / RT = ')
        disp(uu)
        disp('S-S* / R = ')
        disp(ss)
        disp('Θ= ')
        disp(theta)
        disp('φ= ')
        disp(phi)
        
    case 3
        alpha = (1+(0.48+1.574*w-0.176*w^2)*(1-sqrt(Tr)))^2
        A = 0.42748*(Pr/Tr^(2.5))*alpha
        B = 0.08664*(Pr/Tr)
        p = -1
        q = A-B-B^2
        r = -A*B
        z = roots([1 p q r])
        z = z(imag(z)==0)
        
        uu = -3*A/(2*B)*log(1+B/z);
        hh = z-1+uu;
        ss = log(z-B)-A/(2*B)*log(1+B/z);
        theta = (A/B)*log(1+B/z);
        phi = exp(z-1-log(z-B)-theta)
        
        disp('H-H* / RT = ')
        disp(hh)
        disp('U-U* / RT = ')
        disp(uu)
        disp('S-S* / R = ')
        disp(ss)
        disp('Θ= ')
        disp(theta)
        disp('Θ= ')
        disp(theta)
        disp('φ= ')
        disp(phi)
        
    case 4
        alpha = (1+(0.37464+1.54226*w-0.26992*w^2)*(1-sqrt(Tr)))^2
        A = 0.45724*(Pr/Tr^2)*alpha
        B = 0.07780*(Pr/Tr)
        p = -1+B
        q = A-2*B-3*B^2
        r = -A*B+B^2+B^3
        z = roots([1 p q r]);
        z = z(imag(z)==0)
        gamma = (0.37464+1.54226*w-0.26992*w^2)*(sqrt(Tr/alpha));
        
        uu = -A*(1+gamma)/(B*sqrt(8))*log((z+B*(1+sqrt(2)))/(z+B*(1-sqrt(2))));
        hh = z-1+uu;
        ss = log(z-B)-A*gamma/(B*sqrt(8))*log((z+B*(1+sqrt(2)))/(z+B*(1-sqrt(2))));
        theta = (A/(B*sqrt(8))*log((z+B*(1+sqrt(2)))/(z+B*(1-sqrt(2)))));
        phi = exp(z-1-log(z-B)-theta);
        
        disp('H-H* / RT = ')
        disp(hh)
        disp('U-U* / RT = ')
        disp(uu)
        disp('S-S* / R = ')
        disp(ss)
        disp('Θ= ')
        disp(theta)
        disp('φ= ')
        disp(phi)
        
    case 5
        disp('Todo')
        w = input('Introduzca factor acéntrico: \n');
        B0 = 0.083-0.0422/(Tr)^(1.6);
        B1 = 0.139-0.172/(Tr)^(4.2);
        phi = exp(Pr/Tr*(B0+w*B1));
        disp('φ= ')
        disp(phi)
    otherwise
        disp('Te equivocaste chavo\n')
 
end