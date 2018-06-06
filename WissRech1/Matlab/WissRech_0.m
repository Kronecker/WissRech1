n=1000;
h=1/(n+1);


x=[0:n-1]*h;  %f[i]=i*h*(i*h-1);
f=x.*(x-1);

A=eye(n)*2-diag(ones([1,n-1]),1)-diag(ones([1,n-1]),-1);

u=f/A*h*h;

uJacobi=ones(size(u));

lowerDiag=-1;
mainDiag=2;
upperDiag=-1;

fhSquare=f*h*h;
maxIter=100000000;
iteration=0;
lastIterSol=-ones(size(f));
lastIterSol=fhSquare;
actualIteration=lastIterSol*0;
actualIteration=lastIterSol;
while(iteration<maxIter)
    lastIterSol=actualIteration;
    actualIteration(1)=1/mainDiag*(-upperDiag*lastIterSol(2)+fhSquare(1));
   
    actualIteration(2:end-1)=1/mainDiag*(-upperDiag*lastIterSol(3:end)-lowerDiag*lastIterSol(1:end-2)+fhSquare(2:end-1));
   
    actualIteration(end)=1/mainDiag*(-lowerDiag*lastIterSol(end-1)+fhSquare(end));
    
  
   
    if(~mod(iteration,100000))
       % iteration
        norm(abs(lastIterSol)-abs(actualIteration))
        plot(x,u,x,actualIteration);
        drawnow;
    end
    
    iteration=iteration+1;
    
end 
    
    plot(x,u,x,actualIteration);
    
    
    
    
    
    
    
    
