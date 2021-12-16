function z = f(x, Uth)
% x is voicing indicator vector
% Uth : threshold for identify voiced and unvoiced frame
    z = x;
    z(x<Uth)=0;
    
    z(x>1)= 1;
end

