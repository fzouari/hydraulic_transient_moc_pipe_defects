

function [HHm,HHp,QQm,QQp]=solve_transient_moc_pipe_with_defects(dt,H0m,H0p,Q0m,Q0p,HQ1t,HQNt,AAp,aap,kkdb,AAel,ffdp,Nmeas,sys, J, To, EEp, c, orifice)

% transient_moc_sp computes the transient flow (head and flow) in a
% single pipe
%
% Routine written by Fedi Zouari, Hong Kong University of Science and Technology.
% Last modified January 2019
%

global g rho
gamma=rho*g;

if isempty(kkdb); kkdb=zeros(size(AAp)); end
if isempty(AAel); AAel=zeros(size(AAp)); end
if isempty(ffdp); ffdp=zeros(size(AAp)); end

if isvector(H0m)*isvector(H0p)*isvector(Q0m)*isvector(Q0p)*isvector(HQ1t)...
        *isvector(HQNt)*isvector(AAp)*isvector(aap)*isvector(kkdb)*isvector(AAel)~=1
    error('The variables {H0m,H0p,Q0m,Q0p,HQ1t,HQNt,AAp,aap,kkdb,AAel} must be vectors')
end

[H0m,H0p,Q0m,Q0p,AAp,aap,kkdb,AAel,ffdp]= check_size('VEC_ROW',H0m,H0p,Q0m,Q0p,AAp,aap,kkdb,AAel,ffdp);
[HQ1t,HQNt]                             = check_size('VEC_COL',HQ1t,HQNt);

if (length(aap)==1); aap=aap*ones(size(AAp)); end
if (length(ffdp)==1); ffdp=ffdp*ones(size(AAp)); end

check_size('SAME_SIZE',H0m,H0p,Q0m,Q0p,AAp,aap,kkdb,AAel,ffdp);



if nargin>14 
    [J,To] = check_size('VEC_COL',J,To); 
else 
    J=0; EEp=1; c=1; 
    if (length(EEp)==1); EEp=EEp*ones(size(AAp)); end
end

if nargin<19; orifice='ORIF';  end

Nt=length(max(HQNt,HQ1t)); Nx=length(H0m);
              
if length(HQNt)~=Nt; HQNt=HQNt*ones(Nt,1); end
if length(HQ1t)~=Nt; HQ1t=HQ1t*ones(Nt,1); end



AA0  = AAp(1);     AAm  =[AA0,  AAp(1:end-1)];
aa0  = aap(1);     aam  =[aa0,  aap(1:end-1)];
EE0  = EEp(1);     EEm  =[EE0,  EEp(1:end-1)]; 
ffd0 = ffdp(1);    ffdm =[ffd0, ffdp(1:end-1)]; 
dxm  = dt*aam;
dxp  = dt*aap;

QQp = NaN*ones(Nt,length(Nmeas));  % flow rate everywhere at everytime
HHp = NaN*ones(Nt,length(Nmeas));  % head everywhere at everytime
QQm = NaN*ones(Nt,length(Nmeas));  % flow rate everywhere at everytime
HHm = NaN*ones(Nt,length(Nmeas));  % head everywhere at everytime

DDm = sqrt(4*AAm/pi);
DDp = sqrt(4*AAp/pi);
ZZp = aap./(g*AAp);         % Impedance
ZZm = aam./(g*AAm);         % Impedance
RRp = ffdp.*dxp./(2*g*DDp.*AAp.^2); 
RRm = ffdm.*dxm./(2*g*DDm.*AAm.^2); 
DEm = c*gamma*(DDm./EEm);
DEp = c*gamma*(DDp./EEp);

% Initialization of visco-elastic parameters
   CVEp  = zeros(size(AAp)); 
   KVEpj = zeros(size(AAp)); 
   ERKpj = zeros(length(J),length(AAp)); 
   CVEm  = zeros(size(AAp)); 
   KVEmj = zeros(size(AAp)); 
   ERKmj = zeros(length(J),length(AAp)); 


Hdb0=kkdb.*Q0m.^2./(2.*g*AAm.^2); 
% Zdb0=abs(-2*Hdb0./Q1p);

kappa=kkdb./(2*g*AAm.^2);
Cv=1./(2*kappa);
% Cv=(Q1p).^2./(H1p-H1m)/2;

alpha=AAel*sqrt(2*g);
z=zeros(size(AAel)); % elevation of nodes 


% At time j-1: H0 is at time 0, H10 is at time j-1, H1 is at time j
H1p=H0p; H1m=H0m; Q1p=Q0p; Q1m=Q0m; 
H10p=H1p; H10m=H1m; 

% Save Data
QQp(1,:)=Q0p(Nmeas); HHp(1,:)=H0p(Nmeas); 
QQm(1,:)=Q0m(Nmeas); HHm(1,:)=H0m(Nmeas); 

H2p=NaN*ones(1,Nx); Q2p=NaN*ones(1,Nx);
H2m=NaN*ones(1,Nx); Q2m=NaN*ones(1,Nx);

t_start=cputime;
%====================================%
%       Principal Computation        %
%====================================%
for j=2:Nt

   if nargin>15
       CVEp=(aap.^2/g).*DEp.*(sum( J.*(1-exp(-dt./To)) ) ); 
       KVEpj=(aap.^2/g).*DEp.*( (sum( J.*( (dt./To+1).*exp(-dt./To) -1 ) )).*H1p ...
                                -(sum( J.*( (dt./To).*exp(-dt./To) ) )).*H0p ...
                                - 2*(1./DEp).*dt.*(sum( (exp(-dt./To)).*ERKpj./To ,1 )) ); 
       ERKpj=(DEp./2).*( J.*( 1-(To./dt).*(1-exp(-dt./To)) ).*H1p  ...
                         +J.*( -exp(-dt./To) + (To./dt).*(1-exp(-dt./To)) ).*H10p ...
                         +J.*(exp(-dt./To)-1).*H0p ) ...
              +exp(-dt./To).*ERKpj; 
       CVEm=(aam.^2/g).*DEm.*(sum( J.*(1-exp(-dt./To)) ) ); 
       KVEmj=(aam.^2/g).*DEm.*( (sum( J.*( (dt./To+1).*exp(-dt./To) -1 ) )).*H1m ...
                                -(sum( J.*( (dt./To).*exp(-dt./To) ) )).*H0m ...
                                - 2*(1./DEm).*dt.*(sum( (exp(-dt./To)).*ERKmj./To ,1 )) ); 
       ERKmj=(DEm./2).*( J.*( 1-(To./dt).*(1-exp(-dt./To)) ).*H1m  ...
                         +J.*( -exp(-dt./To) + (To./dt).*(1-exp(-dt./To)) ).*H10m ...
                         +J.*(exp(-dt./To)-1).*H0m ) ...
              +exp(-dt./To).*ERKmj; 
    else
       1; 
    end
    
    Cm(1:Nx-1)= ( H1m(2:Nx)-ZZm(2:Nx).*Q1m(2:Nx) - KVEpj(1:Nx-1) )./( 1+CVEp(1:Nx-1) );
    Bm(1:Nx-1)= ( ZZm(2:Nx)+RRm(2:Nx).*abs(Q1m(2:Nx)) )./( 1+CVEp(1:Nx-1) );
    Cp(2:Nx)  = ( H1p(1:Nx-1)+ZZp(1:Nx-1).*Q1p(1:Nx-1) - KVEmj(2:Nx) )./( 1+CVEm(2:Nx) );
    Bp(2:Nx)  = ( ZZp(1:Nx-1)+RRp(1:Nx-1).*abs(Q1p(1:Nx-1)) )./( 1+CVEm(2:Nx) );
    
    %******** Boundary condition i=1
    switch sys
        case 'QH'
        %---- Valve (control of Q)
        Q2p(1)=HQ1t(j); 
        H2p(1)=Cm(1)+Bm(1)*Q2p(1); 
        H2m(1)=H2p(1); Q2m(1)=Q2p(1);
        case 'HQ'
        %---- Reservoir (control of H)
        H2p(1)=HQ1t(j);  
        Q2p(1)=(H2p(1)-Cm(1))/Bm(1); 
        H2m(1)=H2p(1); Q2m(1)=Q2p(1);
    end

    %******** Interior nodes i= 2 : Nx-1
    for i=2:Nx-1
        
        if kappa(i)==0 % ismember(i,Nb)
            Hdb=0;
        else 
            switch orifice
                case 'ORIF'   % Orifice 
            Hdb=(Cm(i)-Cp(i))+ (Bp(i)+Bm(i))^2/(2*kappa(i))* (1 - sqrt( 1-4*kappa(i)*(Cp(i)-Cm(i))/(Bp(i)+Bm(i))^2 ) );
                case 'LIN'    % Linearized 
            % Linearized 1 
            Zdb0=2*kappa.*abs(Q0m); 
            Hdb=-Zdb0(i)*( (Cp-Cm)-(Bm+Bp)*Q0m(i)/2 )/( Bm + Bp + Zdb0(i) ); 
%             % Linearized 2 
%             Hdb=-kappa(i)*abs(Q1p(i))*( Cp-Cm )/( Bm + Bp + kappa(i)*abs(Q1p(i)) ); 
            end
        end
        
        if alpha(i)==0 % ismember(i,Nl)
            QL=0;
        else
            switch orifice
                case 'ORIF'   % Orifice 
            % Orifice 
            QL=-( ( (alpha(i))^2*Bm(i)*Bp(i) )/(2*(Bm(i)+Bp(i))) )* ...
                (1-sqrt(1+(2*(Bm(i)+Bp(i))/(alpha(i)*Bm(i)*Bp(i)))^2*( (Bp(i)*Cm(i)+Bm(i)*Cp(i))/(Bm(i)+Bp(i)) -z(i) ) ) );
                case 'LIN'    % Linearized 
            % Linearized 1 
            QL=alpha(i)*( Cp*Bm + Cm*Bp + (Bm+Bp)*(H0m(i)-2*z(i)) )/( alpha(i)*Bm*Bp + 2*(Bp + Bm)*sqrt(H0m(i)-2*z(i)) ) ; % ( ( ()^2*Bm*Bp )/(2*) )* (1-sqrt(1+(2*(Bm+Bp)/(alpha(i)*Bm*Bp))^2*( (Bp*Cm+Bm*Cp)/(Bm+Bp) -z(i) ) ) );
%             % Linearized 2
%             Zlj=sqrt(H1p-z)./alpha;
%             QL=(1/Zlj(i))*( Cp*Bm + Cm*Bp - (Bm+Bp)*z(i) )/( Bp + Bm + Bm*Bp/Zlj(i) ) ; % ( ( ()^2*Bm*Bp )/(2*) )* (1-sqrt(1+(2*(Bm+Bp)/(alpha(i)*Bm*Bp))^2*( (Bp*Cm+Bm*Cp)/(Bm+Bp) -z(i) ) ) );
%             % Linearized 3 
%             QL=alpha(i)*sqrt( 2*(Cp*Bm + Cm*Bp )/( Bp + Bm) - z(i) ) ; % ( ( ()^2*Bm*Bp )/(2*) )* (1-sqrt(1+(2*(Bm+Bp)/(alpha(i)*Bm*Bp))^2*( (Bp*Cm+Bm*Cp)/(Bm+Bp) -z(i) ) ) );
            end
        end
        
        H2p(i)=(Cp(i)*Bm(i)+Cm(i)*Bp(i))/(Bp(i)+Bm(i))    + (Bm(i)/(Bp(i)+Bm(i)))*Hdb      - (Bp(i)*Bm(i))*QL/(Bp(i)+Bm(i)) ;
        H2m(i)=(Cp(i)*Bm(i)+Cm(i)*Bp(i))/(Bp(i)+Bm(i))    - (Bp(i)/(Bp(i)+Bm(i)))*Hdb      - (Bp(i)*Bm(i))*QL/(Bp(i)+Bm(i)) ;
        Q2p(i)=(Cp(i)-Cm(i))/(Bp(i)+Bm(i))          + (1/(Bp(i)+Bm(i)))*Hdb       - (Bp(i))*QL/(Bp(i)+Bm(i)) ;
        Q2m(i)=(Cp(i)-Cm(i))/(Bp(i)+Bm(i))          + (1/(Bp(i)+Bm(i)))*Hdb       + (Bm(i))*QL/(Bp(i)+Bm(i)) ; 

    end


    %******** Boundary condition i=Nx
    switch sys
        case 'QH'
        %----BC at reservoir
        H2m(Nx)=HQNt(j); 
        Q2m(Nx)=(Cp(Nx)-H2m(Nx))/Bp(Nx);
        H2p(Nx)=H2m(Nx); Q2p(Nx)=Q2m(Nx);
        case 'HQ'
        %----BC at reservoir / upstream valve
        Q2m(Nx)=HQNt(j); 
        H2m(Nx)=Cp(Nx)-Bp(Nx)*Q2m(Nx); 
        H2p(Nx)=H2m(Nx); Q2p(Nx)=Q2m(Nx);
    end
    
        %---- Update and Save
        H10m=H1m; H10p=H1p; % Update Data of time j-1
        Q1m=Q2m; Q1p=Q2p; H1m=H2m; H1p=H2p; % Update Data of time j
        QQm(j,:)=Q1m(Nmeas); QQp(j,:)=Q1p(Nmeas); HHm(j,:)=H1m(Nmeas);  HHp(j,:)=H1p(Nmeas); % Save Data


    clc
    fprintf('\t MOC: \n\n')    
    show_clock(j/Nt,t_start,'m');

end

end




