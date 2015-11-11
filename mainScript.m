clear,close;

%% Matrix Intializations.
% The size of the block
N=8;                        
M=8;
%Read Images
I=imread('cameraman.tif');             
[h,w]=size(I);                      
Blocks(h/N,w/M).block=zeros(N,M);  

Norm_Mat=[16 11 10 16 24 40 51 61       
          12 12 14 19 26 58 60 55
          14 13 16 24 40 57 69 56
          14 17 22 29 51 87 80 62
          18 22 37 56 68 109 103 77
          24 35 55 64 81 104 113 92
          49 64 78 87 103 121 120 101
          72 92 95 98 112 100 103 99];
      
%% DCT Part
for a = 1:h/N
    for b = 1:w/M
        for u = 1:N
            for v = 1:M
                sum = 0;
                for i=1:N
                    for j=1:M
                        sum = sum+ double(I(N*(a-1)+i,M*(b-1)+j))*cos(pi*(2*i+1)*u/(2*N))*cos(pi*(2*j+1)*v/(2*M));
                    end
                end
                if u==0
                    sum=sum*(1/sqrt(2));
                else
                    sum=sum*1;
                end
                if v==0
                    sum=sum*(1/sqrt(2));
                else
                    sum=sum*1;
                end
                Blocks(a,b).block(u,v) = sum*1/4;
            end
        end
        %Blocks(a,b).block=round(Blocks(a,b).block./Norm_Mat);
    end
end

%% Inverse DCT Part
DecodeImg = zeros(w,h);
for a = 1:h/N
    for b = 1:w/M
        for u = 1:N
            for v = 1:M
                sum = 0;
                for i=1:N
                    for j=1:M
                        tmpSum = double(Blocks(a,b).block(i,j))*cos(pi*(2*i+1)*u/(2*N))*cos(pi*(2*j+1)*v/(2*M));
                        if u==0
                            tmpSum=tmpSum*(1/sqrt(2));
                        else
                            tmpSum=tmpSum*1;
                        end
                        if v==0
                            tmpSum=tmpSum*(1/sqrt(2));
                        else
                            tmpSum=tmpSum*1;
                        end
                        sum = sum+tmpSum;
                    end
                end
                DecodeImg(N*(a-1)+u,M*(b-1)+v) = sum*(1/4);
            end
        end        
    end
end
imshow(DecodeImg);