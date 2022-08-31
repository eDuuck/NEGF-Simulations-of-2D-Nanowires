clear all
hbar =1.06e-34; q =1.6e-19; m =0.1*9.1e-31; qh=q/ hbar ;
a =2.5e-9; t0 =( hbar ^2) /(2* m*(a^2)*q);
NW =25; Np =5; L= zeros(Np);R=L;L(1 ,1)=1; R(Np ,Np)=1;
zplus =i*1e-12;

% Hamiltonian

al =4* t0;by=-t0;bx=-t0;
alpha = kron (eye(NW),al)+ kron(diag(ones(1,NW - 1) ...
,+1) ,by)+ kron ( diag ( ones (1,NW -1) ,-1),by');
alpha = alpha + diag ([1:1: NW ].*0) ;

EE=t0;ii =0; for B =20
    ii=ii +1; E(ii)=B;
    ig0 =( EE+ zplus )*eye(NW)-alpha ;
    if ii ==1
        gs1=inv(ig0);gs2=inv(ig0);end

    beta = kron ( diag ( exp (i*qh*B*a*a *[1:1: NW ])),bx);
    H= kron (eye(Np),alpha );
    if Np >1
        H=H+ kron ( diag ( ones (1,Np -1) ,+1) ,beta )+ kron ( diag ( ones ....
            (1,Np -1) ,-1),beta'); end

    change = 1;
    while change >5e-5
        Gs=inv (ig0 -beta'* gs1 * beta);
        change =sum (sum (abs (Gs -gs1 )))/( sum (sum (abs (gs1 )+abs(...
        Gs))));
        gs1 =0.5* Gs +0.5* gs1;
    end 
    sig1 =beta'* gs1 * beta ; sig1 = kron (L, sig1 ); gam1 =i*( sig1 -sig1');
    %subplot(1,3,1);imagesc(abs(sig1));subplot(1,3,2);imagesc(real(sig1));subplot(1,3,3);imagesc(imag(sig1));
    %pause(0.1);

    change = 1;
    while change >5e-5
        Gs=inv (ig0 -beta* gs2 * beta');
        change =sum (sum (abs (Gs -gs2 )))/( sum (sum (abs (gs2 )+abs(...
        Gs))));
        gs2 =0.5* Gs +0.5* gs2;
    end
    sig2 =beta* gs2 * beta'; sig2 = kron (R, sig2 ); gam2 =i*( sig2 -sig2');

    G=inv (( EE*eye(Np*NW))-H-sig1 - sig2 );
    Gn=G* gam1 *G';


    A=i*(G-G'); V= real ( diag (Gn ./A));
    Tcoh = real ( trace ( gam1 *G* gam2 *G'));TM=real(trace ( gam2*Gn));
    imagesc(V);
    pause(0.1);
    Y(ii) = (V(1) - V(NW))/Tcoh;ii
end

hold on
h= plot (E,Y,'k');
set(h,'linewidth' ,[3.0])
set(gca ,'Fontsize' ,[36])
xlabel ('B- field (T) --->')
ylabel ('R_{xy} --->')
grid on