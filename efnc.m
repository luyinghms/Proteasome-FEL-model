function ret=efnc(x,para,flag)
%evaluation function: 1: gaussian function; 2: exponential function. 3:
%step function

switch flag
    case 1
        ret=exp(-(x-para(1))^2/2/para(2)^2);
    case 2
       if x>0
            ret=exp(-x/para(1));
       else 
            ret=1;
       end
          
    case 3
        ret=exp(para(3)*(x-para(1))/para(2))/(1+exp(para(3)*(x-para(1))/para(2)));
end
        
end
