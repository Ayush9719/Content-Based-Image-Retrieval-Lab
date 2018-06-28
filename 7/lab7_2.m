fr=zeros(100,1);
srcFiles = dir('*.jpg');
Fdirectiond = zeros(100,1);
Fcontrastd = zeros(100,1);
Flind = zeros(100,1);
Fregd= zeros(100,1);
Fcoarsenessd= zeros(100,1);
Froughnessd = zeros(100,1);
name = zeros(100,1);
for i=1:100
    name(i) = i;
    d = strcat(srcFiles(i).name);
    IColor = imread(d);
    I = rgb2gray(IColor);
    [r,c] = size(I);
    G=double(I);
    [counts,graylevels]=imhist(I);
    PI=counts/(r*c);
    averagevalue=sum(graylevels.*PI);
    u4=sum((graylevels-repmat(averagevalue,[256,1])).^4.*PI);
    variance=sum((graylevels-repmat(averagevalue,[256,1])).^2.*PI);
    alpha4=u4/variance^2;
    Fcontrastd(i)=sqrt(variance)/alpha4.^(1/4);



    A1=zeros(r,c);A2=zeros(r,c);
    A3=zeros(r,c);A4=zeros(r,c);
    A5=zeros(r,c);A6=zeros(r,c);
    Sbest=zeros(r,c);

    E1h=zeros(r,c);E1v=zeros(r,c);
    E2h=zeros(r,c);E2v=zeros(r,c);
    E3h=zeros(r,c);E3v=zeros(r,c);
    E4h=zeros(r,c);E4v=zeros(r,c);
    E5h=zeros(r,c);E5v=zeros(r,c);
    E6h=zeros(r,c); E6v=zeros(r,c);
    flag=0;    

    for x=2:r
    for y=2:c
                A1(x,y)=(sum(sum(G(x-1:x,y-1:y))));
    end
    end
    for x=2:r-1
    for y=2:c-1
                E1h(x,y) = A1(x+1,y)-A1(x-1,y);
                E1v(x,y) = A1(x,y+1)-A1(x,y-1);
    end
    end
        E1h=E1h/2^(2*1);
        E1v=E1v/2^(2*1);


    if (r<4||c<4)
            flag=1;
    end

    if(flag==0)
    for x=3:r-1
    for y=3:c-1
                    A2(x,y)=(sum(sum(G(x-2:x+1,y-2:y+1))));
    end
    end
    for x=3:r-2
    for y=3:c-2
                    E2h(x,y) = A2(x+2,y)-A2(x-2,y);
                    E2v(x,y) = A2(x,y+2)-A2(x,y-2);
    end
    end
    end
        E2h=E2h/2^(2*2);
        E2v=E2v/2^(2*2);

    if (r<8||c<8)
            flag=1;
    end

    if(flag==0)
    for x=5:r-3
    for y=5:c-3
                    A3(x,y)=(sum(sum(G(x-4:x+3,y-4:y+3))));
    end
    end
    for x=5:r-4
    for y=5:c-4
                    E3h(x,y) = A3(x+4,y)-A3(x-4,y);
                    E3v(x,y) = A3(x,y+4)-A3(x,y-4);
    end
    end
    end
        E3h=E3h/2^(2*3);
        E3v=E3v/2^(2*3);


    if (r<16||c<16)
            flag=1;
    end

    if(flag==0)
    for x=9:r-7
    for y=9:c-7
                    A4(x,y)=(sum(sum(G(x-8:x+7,y-8:y+7))));
    end
    end
    for x=9:r-8
    for y=9:c-8
                    E4h(x,y) = A4(x+8,y)-A4(x-8,y);
                    E4v(x,y) = A4(x,y+8)-A4(x,y-8);
    end
    end
    end
        E4h=E4h/2^(2*4);
        E4v=E4v/2^(2*4);


    if (r<32||c<32)
            flag=1;
    end

    if(flag==0)
    for x=17:r-15
    for y=17:c-15
                    A5(x,y)=(sum(sum(G(x-16:x+15,y-16:y+15))));
    end
    end
    for x=17:r-16
    for y=17:c-16
                    E5h(x,y) = A5(x+16,y)-A5(x-16,y);
                    E5v(x,y) = A5(x,y+16)-A5(x,y-16);
    end
    end
    end
        E5h=E5h/2^(2*5);
        E5v=E5v/2^(2*5);


    if (r<64||c<64)
            flag=1;
    end

    if(flag==0)
    for x=33:r-31
    for y=33:c-31
                    A6(x,y)=(sum(sum(G(x-32:x+31,y-32:y+31))));
    end
    end
    for x=33:r-32
    for y=33:c-32
                    E6h(x,y) = A6(x+32,y)-A6(x-32,y);
                    E6v(x,y) = A6(x,y+32)-A6(x,y-32);
    end
    end
    end
        E6h=E6h/2^(2*6);
        E6v=E6v/2^(2*6);

    for ii=1:r
        for j=1:c
                    [maxv,index]=max([abs(E1h(ii,j)),abs(E1v(ii,j)),abs(E2h(ii,j)),abs(E2v(ii,j)),...
                        abs(E3h(ii,j)),abs(E3v(ii,j)),abs(E4h(ii,j)),abs(E4v(ii,j)),abs(E5h(ii,j)),...
                        abs(E5v(ii,j)),abs(E6h(ii,j)),abs(E6v(ii,j))]);
                    k=floor((index+1)/2);

                    Sbest(ii,j)=2.^k;
        end
    end

    Fcoarsenessd(i)=sum(sum(Sbest))/(r*c);

    PrewittH = [-1 0 1;-1 0 1;-1 0 1];
    PrewittV = [1 1 1;0 0 0;-1 -1 -1];

    deltaH=zeros(r,c);
    for ii=2:r-1
        for j=2:c-1
                    deltaH(ii,j)=sum(sum(G(ii-1:ii+1,j-1:j+1).*PrewittH));
        end
    end
    for j=2:c-1
            deltaH(1,j)=G(1,j+1)-G(1,j);
            deltaH(r,j)=G(r,j+1)-G(r,j);
    end
    for ii=1:r
            deltaH(ii,1)=G(ii,2)-G(ii,1);
            deltaH(ii,c)=G(ii,c)-G(ii,c-1);
    end
        deltaV=zeros(r,c);
    for ii=2:r-1
        for j=2:c-1
                    deltaV(ii,j)=sum(sum(G(ii-1:ii+1,j-1:j+1).*PrewittV));
        end
    end
    for j=1:c
            deltaV(1,j)=G(2,j)-G(1,j);
            deltaV(r,j)=G(r,j)-G(r-1,j);
    end
    for ii=2:r-1
            deltaV(ii,1)=G(ii+1,1)-G(ii,1);
            deltaV(ii,c)=G(ii+1,c)-G(ii,c);
    end
        deltaG=(abs(deltaH)+abs(deltaV))/2;
        theta=zeros(r,c);
    for ii=1:r
        for j=1:c
            if (deltaH(ii,j)==0)&&(deltaV(ii,j)==0)
                            theta(ii,j)=0;
            elseif deltaH(ii,j)==0
                            theta(ii,j)=pi;
            else
                            theta(ii,j)=atan(deltaV(ii,j)/deltaH(ii,j))+pi/2;
            end
        end
    end
    deltaGt = deltaG(:);
    theta1=theta(:);

    n = 16;
    HD = zeros(1,n);
    Threshold=12;
    counti=0;
    for m=0:(n-1)
        countk=0;
        for k = 1:length(deltaGt)
            if ((deltaGt(k)>=Threshold) && (theta1(k)>=(2*m-1)*pi/(2*n)) && (theta1(k)<(2*m+1)*pi/(2*n)))
                        countk=countk+1;
                        counti=counti+1;
            end
        end
        HD(m+1) = countk;
    end
    HDf = HD/counti;
    [m, p]=findpeaks(HDf,0.000005);

    Fd=0;
    for np = 1:length(m)
            phaiP=m(np)*(pi/n);
        for phi=1:length(HDf)
                    Fd=Fd+(phi*(pi/n)-phaiP)^2*HDf(phi);
        end
    end
    lamda = 0.03;
    Fdirectiond(i) = 1 - lamda*np*Fd;

    Froughnessd(i) = Fcontrastd(i) + Fcoarsenessd(i);
    
    glcm3 = graycomatrix(I,'Offset',[1 1]);
    num = 0;
    den = 0;
    for a=1:8
        for b=1:8
            num = num + glcm3(a,b)*cos((a-b)*2*pi/8);
            den = den + glcm3(a,b);
        end
    end
    Flind(i) = num/den;
    Fregd(i) = 1 -lamda*(Fcoarsenessd(i) + Fcontrastd(i) + Fdirectiond(i) + Flind(i));
end
res = [name,abs(Fcoarsenessd),abs(Fcontrastd),abs(Fdirectiond),abs(Flind),abs(Fregd),abs(Froughnessd)];
head = {'FileName','Coarsness','Contrast','Directionality','Line-Likeliness','Regularity','Roughness'};
xlswrite('tamura',head,'Sheet1');
xlswrite('tamura',res,'Sheet1');