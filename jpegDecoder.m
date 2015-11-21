%% Decoding Part
clear;
N=8;
M=8;
Norm_Mat=[16 11 10 16 24 40 51 61       
          12 12 14 19 26 58 60 55
          14 13 16 24 40 57 69 56
          14 17 22 29 51 87 80 62
          18 22 37 56 68 109 103 77
          24 35 55 64 81 104 113 92
          49 64 78 87 103 121 120 101
          72 92 95 98 112 100 103 99];
%% Import TXT data
fileID = fopen('code.txt');
%load('matlab.mat');
C = textscan(fileID,'%d');
code = C{1}';
fclose(fileID);

%% Parse information
h = code(1); w = code(2);
padH= code(3); padW = code(4);
DCcoes = code(5:4+(h/N)*(w/M));
RunLengthArray = double(code(5+(h/N)*(w/M):end));
DeBlocks(h/N,w/M).block=zeros(N,M);
ZigZagBlocks(h/N,w/M).array=zeros(1,N*M);

%% Run Length Decoding
RunLengthIndex = 1;
for a=1:h/N
    for b=1:w/M
        ZigZagBlocks(a,b).array=zeros(1,N*M);
        ZigZagBlockIndex = 1;
        
        while RunLengthArray(RunLengthIndex+1) ~= 0  
            ZigZagBlockIndex = ZigZagBlockIndex + RunLengthArray(RunLengthIndex)+1;
            ZigZagBlocks(a,b).array(ZigZagBlockIndex) = double(RunLengthArray(RunLengthIndex+1));
            RunLengthIndex = RunLengthIndex+2;
            if ZigZagBlockIndex == N*M
                break
            end
        end
        if RunLengthArray(RunLengthIndex+1) == 0
            RunLengthIndex = RunLengthIndex+2;
        end
    end
end
%% Inverse DC and ZigZag reconstruction calculations
DCindex = 1;
PreDC = 0;
for a=1:h/N
    for b=1:w/M
        ZigZagBlocks(a,b).array(1) = double(DCcoes(DCindex)+PreDC);
        PreDC = ZigZagBlocks(a,b).array(1);
        DCindex = DCindex + 1;
        %%ZigZag reconstruction
        DeBlock=zeros(N,N);
        count=1;
        for s=1:N
            if mod(s,2)~=0
                for m=s:-1:1
                    DeBlock(m,s+1-m)= ZigZagBlocks(a,b).array(count);
                    count=count+1;
                end;
            else
                for m=1:s
                    DeBlock(m,s+1-m)= ZigZagBlocks(a,b).array(count);
                    count=count+1;
                end;
            end;
        end;
        if mod(N,2)~=0
            flip=1;
        else
            flip=0;
        end;
        for s=N+1:2*N-1
            if mod(flip,2)~=0
                for m=N:-1:s+1-N
                    DeBlock(m,s+1-m)= ZigZagBlocks(a,b).array(count);
                    count=count+1;
                end;
            else
                for m=N:-1:s+1-N
                    DeBlock(s+1-m,m)= ZigZagBlocks(a,b).array(count);
                    count=count+1;
                end;
            end;
            flip=flip+1;
        end;
        DeBlocks(a,b).block = double(double(DeBlock.*Norm_Mat));
    end
end

%% Inverse DCT Part
DecodeImg = zeros(h,w);
for a = 1:h/N
    for b = 1:w/M
        for u = 1:N
            for v = 1:M
                sum = 0;
                for i=1:N
                    for j=1:M
                        tmpSum = double(DeBlocks(a,b).block(i,j))*...
                            cos(pi*(2*(u-1)+1)*(i-1)/(2*N))*...
                            cos(pi*(2*(v-1)+1)*(j-1)/(2*M));
                        
                        if (i-1)==0
                            tmpSum=tmpSum*sqrt(1/2);
                        else
                            tmpSum=tmpSum*1;
                        end
                        if (j-1)==0
                            tmpSum=tmpSum*sqrt(1/2);
                        else
                            tmpSum=tmpSum*1;
                        end
                        sum = sum+tmpSum;
                    end
                end
                DecodeImg(N*(a-1)+u,M*(b-1)+v) = double(sum*(1/4));
            end
        end        
    end
end
finalImage = DecodeImg(1:h-padH,1:w-padW);
imshow(uint8(finalImage));