function [samples, n_samples] = load_samples(samplefilepath, format, little_endian, nr)
 %%%%%%%%%%%%%%%%%%%%%%
 %% read the samples
 %%%%%%%%%%%%%%%%%%%%%%

 if nargin < 4
     nr = inf;
 end
 if nargin > 1 && strcmp(format,'cplx')
   fid = fopen(samplefilepath, 'r');
   if (nargin < 3)
       little_endian = false;
   end
   if (little_endian)
       [s, n] = fread(fid, nr, 'int32', 0, 'ieee-le');
   else
       [s, n] = fread(fid, nr, 'int32', 0, 'ieee-be');
   end
   fclose(fid);
   s = reshape(s, 2, []).';
   s = s(:,1) + 1i * s(:,2);
   n = n/2;
 end

 if nargin > 1 && strcmp(format,'float32')
   fid = fopen(samplefilepath, 'r');
   [s, n] = fread(fid, nr, 'float');
   fclose(fid);
   s = reshape(s, 2, []).';
   s = s(:,1) + 1i * s(:,2);
   n = n/2;
 end

 s(1:10,:)
 samples = s;
 n_samples = n;
end