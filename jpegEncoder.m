clear,close;

%% Matrix Intializations.
% The size of the block
N=8;                        
M=8;
%Read Images
I=imread('cameraman.tif'); 
[h,w]=size(I);                      
Blocks(h/N,w/M).block=zeros(N,M);

%%Padding the image if needed
modN = mod(h,N); modM = mod(w,M);
gapN = N;gapM = M;
if modN == 0 
    modN = N; 
    gapN = 0;
end
if modM == 0
    modM = M;
    gapM = 0;
end
newH = h+(N-modN);newW = w+(M-modM);
padH = round((N-modN)/2); 
padW = round((M-modM)/2);
%Make the zero padding symmetrical distributing around the image 
paddingI=[zeros(padH,newW);
    zeros(h,padW),I,zeros(h,gapM-padW);
    zeros(gapN-padH,newW)];
h = newH; w = newW; I = paddingI;

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
                        sum = sum+ double(I(N*(a-1)+i,M*(b-1)+j))*...
                            cos(pi*(2*(i-1)+1)*(u-1)/(2*N))*...
                            cos(pi*(2*(j-1)+1)*(v-1)/(2*M));
                    end
                end
                if (u-1)==0
                    sum=sum*sqrt(1/2);
                else
                    sum=sum*1;
                end
                if (v-1)==0
                    sum=sum*sqrt(1/2);
                else
                    sum=sum*1;
                end
                Blocks(a,b).block(u,v) = sum*(1/4); 
            end
        end
        Blocks(a,b).block=round(Blocks(a,b).block./double(Norm_Mat));
    end
end

%% DC coefficiences calculations && ZigZag Recordering
DCcoes = zeros(1,length(Blocks));
DCindex = 1; PreDC = 0;
DCcoes(1) = Blocks(1,1).block(1,1);
ZigZagBlocks(h/N,w/M).array = zeros(1,N*M);
for a=1:h/N
    for b=1:w/M
        DCcoes(DCindex) = Blocks(a,b).block(1,1) - PreDC;
        PreDC = Blocks(a,b).block(1,1);
        DCindex = DCindex + 1;
        %%Zig Zag Part
        ZigZagArray = zeros(1,N*M);
        count=1;
        for s=1:N
            if mod(s,2)~=0
                for m=s:-1:1
                    ZigZagArray(count)=Blocks(a,b).block(m,s+1-m);
                    count=count+1;
                end;
            else
                for m=1:s
                    ZigZagArray(count)=Blocks(a,b).block(m,s+1-m);
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
                    ZigZagArray(count)=Blocks(a,b).block(m,s+1-m);
                    count=count+1;
                end;
            else
                for m=N:-1:s+1-N
                    ZigZagArray(count)=Blocks(a,b).block(s+1-m,m);
                    count=count+1;
                end;
            end;
            flip=flip+1;
        end;
        ZigZagBlocks(a,b).array = double(ZigZagArray);
    end
end


%% Run Length Codding for ZigZagBlocks
RunLengthArray = zeros(1,0);
for a=1:h/N
    for b=1:w/M
        zeroCount = 0;
        for i = 2:length(ZigZagBlocks(a,b).array)
            if ZigZagBlocks(a,b).array(i) == 0 
                zeroCount = zeroCount+1;
            else
                index = length(RunLengthArray) + 1;
                RunLengthArray(index) = zeroCount;
                RunLengthArray(index+1) = double(ZigZagBlocks(a,b).array(i));
                zeroCount = 0;
            end
        end
        if zeroCount ~=0
            index = length(RunLengthArray) + 1;
            RunLengthArray(index) = 0;
            RunLengthArray(index+1) = 0;
        end
    end
end
% %% Output TXT file with format: [height,width,DCCoes(with length of (height/8 * width/8)),RunLenngth array with format(01,x1,02,x2..)]
code = [h,w,DCcoes,RunLengthArray];
fileID = fopen('code.txt','w');
fprintf(fileID,'%d ',code);
fclose(fileID);

