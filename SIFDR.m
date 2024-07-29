
%¡°A real-time unsupervised hyperspectral band selection via spatial-spectral information fusion based downscaled region¡±

% Input:
% pre     - Raw image data(after bicubic interpolation), size=Np*Nb, note: Np is the number of pixels, Nb is the number of bands.
% bandnum - Number of band requests


% Output:
% C       - Band index

function [C] = SIFDR(pre,bandnum)
[Np,Nb] = size(pre);
%% +++++++++++++++++Spectral Information Estimation+++++++++++++++%%
R=(1/Np)*(pre'*pre);
d=pre';
X=inv(R);
group=ceil(Np/1000);
w=[];
for i = 1:group
up=min((1000*i),Np);
do=1000*(i-1)+1;
di = d(:,do:up); 
wi=(X*di)*inv(di'*X*di);                     % Weight vector of bands from Eq.6.
w=[w,wi];
end
w=abs(w);
rho=mean(w,2);                               % Spectral information estimation value "¦Ñ" of bands from Eq.7.


%% +++++++++++++++++Band Weak Redundancy Ranking+++++++++++++++++ %%
%% delta calculation 
Dist_matrix= (G_D(pre))/sqrt(Np);            % Tanh normalization from Eq.14.
Dist_matrix=abs(tanh((Dist_matrix).^1));
[~,ordrho]=sort(rho,'descend');
delta = zeros(Nb,1);
delta(ordrho(1))=-1;
maxD = max(max(Dist_matrix));
for i = 2:Nb
    delta(ordrho(i))=maxD;
    for j=1:i-1
        if Dist_matrix(ordrho(i),ordrho(j))<delta(ordrho(i))
            delta(ordrho(i)) = Dist_matrix(ordrho(i),ordrho(j));           % "¦Ä" calculation by Eq.10.
        end
    end
end
delta(ordrho(1)) = max(delta);
rho = (rho-min(rho(:)))/(max(rho(:))-min(rho(:)));                         % "¦Ñ" normalization from Eq.13.

%% The final importance
phi =rho.*delta;                          % Rangking metric"¦Õ" calculation by Eq.12.

%% Bands Ranking and selecting
[~,order_band] = sort(phi,'descend');
C = order_band(1:bandnum);
end

