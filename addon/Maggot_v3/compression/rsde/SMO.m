function alpha=SMO(Q,D)
%
% Sequential Minimal Optimisation (SMO) algorithm for Reduced Set Density Estimation (RSDE).
%    Minimising 0.5*alpha*Q*alpha'-alpha*D'
%     
%    Use format: alpha=SMO(Q,D)
%
%    Input:  Q [NxN]:        Kernal matrix 
%            D [1xN]:        Parzen density estimate   
%    Return: alpha  [1xN]:   Weight vector obtained by RSDE
%
%    Technical reference: B. Scholkopf, J. Platt and J. Shawe-Taylor et al. "Estimating the support
%            of a high-dimensional distribution", Neural Computation, 13: 1443-1471, 2001.
%
%    Copyright Chao He & Mark Girolami
%    Last revised on January 22th, 2003
%    Acknowledgement: Thanks to Anna Szymkowiak from Technical University of Denmark 
%                     for vectorising certain parts of the code.


%Initialisation
alpha=D./sum(D);

examineAll=1;
alpha_tolerance= 1e-6;  %Tolerance to threshold weight be zero     
error_tolerance= 1e-8;  %Iteration error tolerance to terminate the algorithm

E2=[];
stop=0;
loop=0;
loopnum=0;
while (~stop)
    if (examineAll)
        loopnum=loopnum+1;
     %   fprintf('Loop %d ...\n',loopnum);
        s=find(alpha>alpha_tolerance);
        alpha_tmp=alpha(s);
        [alpha_max,I_max]=max(alpha_tmp);
        I2=s(I_max);
        numChanged=0;
        [alpha,I1,no]=searchPoint(I2,alpha,Q,D);    
        numChanged=numChanged+no;
        tt=find(alpha_tmp~=alpha_max&alpha_tmp~=alpha(I1));
        s=s(tt);
        if length(s)==0
            loop=1;
        end
        examineAll=0;
    else 
        alpha_tmp=alpha(s);
        [alpha_max,I_max]=max(alpha_tmp);
        I2=s(I_max);
        alpha_bk=alpha;
        [alpha,I1,no]=searchPoint(I2,alpha,Q,D);    
        numChanged=numChanged+no;
        ddd=find(alpha_tmp~=alpha_max&alpha_tmp~=alpha(I1));
        s=s(ddd);
        if length(s)==0
            loop=1;
        end
    end
    if loop==1
        E1=0.5*alpha*Q*alpha'-alpha*D';
        %fprintf('Iteration error: %e\n\n',E1-E2);
        
        E=[E2,E1];
        if E1>E2
            alpha=alpha_bk;     
            stop=1;
        else if (abs(E1-E2)<error_tolerance)
                stop=1;
            end
        end
        examineAll=1;
        loop=0; 
        E2=E1;
        if numChanged==0
            stop=1;
        else
            examineAll=1;
            loop=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,I1,numChanged]=searchPoint(I2,alpha,Q,D)
% Find the second point to be updated
alpha_tolerance=1e-6;
[N,N]=size(Q);

%%%%%%%%%%%%%%%%%%%
dW=zeros(1,N);
W1=alpha*Q-D;
W2=repmat(alpha*Q(I2,:)'-D(I2),1,N);
ind=find(alpha>alpha_tolerance);
dummy=abs(W1-W2);
dW(1,ind)=dummy(ind);
%%%%%%%%%%%%%%%%%%%%


[dW_max,I1]=max(dW);
dW=W1-W2;
[alpha,numChanged]=updateWeight(I1,I2,alpha,dW(I1),Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,numChanged]=updateWeight(I1,I2,alpha,dW,Q)
%Updating the weights
alpha_tolerance=1e-6;
if I1==I2
    numChanged=0;
    return;
end
if dW==0
    numChanged=0;
    return;
end
if alpha(I1)<alpha_tolerance
    numChanged=0;
    return;
end
alph2=alpha(I2)+dW/(Q(I1,I1)-2*Q(I1,I2)+Q(I2,I2));
if alph2<0
    alph2=0;
end
alph1=alpha(I1)+alpha(I2)-alph2;
if alph1<0
    alph1=0;
    alph2=alpha(I1)+alpha(I2);
end
alpha(I1)=alph1;
alpha(I2)=alph2;
numChanged=1; 
%fprintf('I2=%d,  alpha2=%f\n',I2,alpha(I2));
%fprintf('I1=%d,  alpha1=%f\n\n',I1,alpha(I1));

    