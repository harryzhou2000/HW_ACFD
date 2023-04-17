function [rhoL,uL,pL,rhoR,uR,pR] = getRiemannProblem(pnum)



switch pnum
    case 1
        rhoL = 1;       uL = 0.75;      pL = 1   ;    rhoR = 0.125;   uR = 0;        pR = 0.1;    
    case 2
        rhoL = 1;       uL = -2  ;      pL = 0.4 ;    rhoR = 1    ;   uR = 2;        pR = 0.4;  
    case 3
        rhoL = 1;       uL = 0   ;      pL = 1000;    rhoR = 1    ;   uR = 0;        pR = 0.01;   
    case 4
        rhoL = 5.99924; uL = 19.5975;   pL = 460.894; rhoR = 5.99242; uR = -6.19633; pR = 46.0950;
    case 5
        rhoL = 1;       uL = -19.59745; pL = 1000;    rhoR = 1    ;   uR = -19.59745;pR = 0.01;   
    otherwise
        error(' ');
end

  


