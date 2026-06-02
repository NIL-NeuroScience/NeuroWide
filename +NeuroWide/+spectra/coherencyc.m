% CHRONUX coherencyc function
% https://chronux.org
function [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
% Multi-taper coherency,cross-spectrum and individual spectra - continuous process
% 
% Usage:
% [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       C (magnitude of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi (phase of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S12 (cross spectrum -  frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S1 (spectrum 1 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S2 (spectrum 2 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - theoretical/jackknife (depending on err(1)=1/err(1)=2) standard deviation for phi. 
%                Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
%                bands for phi - only for err(1)>=1 
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 2; error('Need data1 and data2'); end;
data1=change_row_to_column(data1);
data2=change_row_to_column(data2);
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave]=getparams(params);
if nargout > 8 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 6 && err(1)==0;
%   Errors computed only if err(1) is nonzero. Need to change params and run again.
    error('When errors are desired, err(1) has to be non-zero.');
end;
N=check_consistency(data1,data2);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J1=mtfftc(data1,tapers,nfft,Fs);
J2=mtfftc(data2,tapers,nfft,Fs);
J1=J1(findx,:,:); J2=J2(findx,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
C=abs(C12); 
phi=angle(C12);
if nargout>=9; 
     [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave);
elseif nargout==8;
     [confC,phistd]=coherr(C,J1,J2,err,trialave);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=change_row_to_column(data)
% Helper routine to transform 1d arrays into column vectors that are needed
% by other routines in Chronux
%
% Usage: data=change_row_to_column(data)
% 
% Inputs:
% data -- required. If data is a matrix, it is assumed that it is of the
% form samples x channels/trials and it is returned without change. If it
% is a vector, it is transformed to a column vector. If it is a struct
% array of dimension 1, it is again returned as a column vector. If it is a
% struct array with multiple dimensions, it is returned without change
% Note that the routine only looks at the first field of a struct array.
% 
% Ouputs:
% data (in the form samples x channels/trials)
%

dtmp=[];
if isstruct(data);
   C=length(data);
   if C==1;
      fnames=fieldnames(data);
      eval(['dtmp=data.' fnames{1} ';'])
      data=dtmp(:);
   end
else
  [N,C]=size(data);
  if N==1 || C==1;
    data=data(:);
  end;
end;
end

function [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
% Helper function to convert structure params to variables used by the
% various routines - also performs checks to ensure that parameters are
% defined; returns default values if they are not defined.
%
% Usage: [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
%
% Inputs:
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%           - optional
%             tapers : precalculated tapers from dpss or in the one of the following
%                       forms:  
%                       (1) A numeric vector [TW K] where TW is the
%                           time-bandwidth product and K is the number of
%                           tapers to be used (less than or equal to
%                           2TW-1). 
%                       (2) A numeric vector [W T p] where W is the
%                           bandwidth, T is the duration of the data and p 
%                           is an integer such that 2TW-p tapers are used. In
%                           this form there is no default i.e. to specify
%                           the bandwidth, you have to specify T and p as
%                           well. Note that the units of W and T have to be
%			                consistent: if W is in Hz, T must be in seconds
% 			                and vice versa. Note that these units must also
%			                be consistent with the units of params.Fs: W can
%		    	            be in Hz if and only if params.Fs is in Hz.
%                           The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Outputs: 
% The fields listed above as well as the struct params. The fields are used
% by some routines and the struct is used by others. Though returning both
% involves overhead, it is a safer, simpler thing to do.

if ~isfield(params,'tapers') || isempty(params.tapers);  %If the tapers don't exist
     display('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers=[3 5];
end;
if ~isempty(params) && length(params.tapers)==3 
    % Compute timebandwidth product
    TW = params.tapers(2)*params.tapers(1);
    % Compute number of tapers
    K  = floor(2*TW - params.tapers(3));
    params.tapers = [TW  K];
end

if ~isfield(params,'pad') || isempty(params.pad);
    params.pad=0;
end;
if ~isfield(params,'Fs') || isempty(params.Fs);
    params.Fs=1;
end;
if ~isfield(params,'fpass') || isempty(params.fpass);
    params.fpass=[0 params.Fs/2];
end;
if ~isfield(params,'err') || isempty(params.err);
    params.err=0;
end;
if ~isfield(params,'trialave') || isempty(params.trialave);
    params.trialave=0;
end;

tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;
err=params.err;
trialave=params.trialave;
end

function [N,C]=check_consistency(data1,data2,sp)
% Helper routine to check consistency of data dimensions
% Usage: [N,C]=check_consistency(data1,data2,sp)
% Inputs:
% data1 - first dataset
% data2 - second dataset
% sp - optional argument to be input as 1 when one of the two data sets is
% spikes times stored as a 1d array.
% Outputs:
% Dimensions of the datasets - data1 or data2 (note that 
%    routine stops with an error message if dimensions don't match - [N,C]
%    N left empty for structure arrays
N1=[]; N2=[];
if nargin < 3 || isempty(sp); sp=0; end;
if isstruct(data1);
    C1=length(data1);
else
    [N1,C1]=size(data1);
end;
if isstruct(data2);
    C2=length(data2);
else
    [N2,C2]=size(data2);
end;
if C1~=C2; error('inconsistent dimensions'); end;
if sp==0;
   if ~isstruct(data1) && ~isstruct(data2);
      if N1~=N2; error('inconsistent dimensions'); end;
   end;
end;
N=N1; C=C1;
end

function [f,findx]=getfgrid(Fs,nfft,fpass)
% Helper function that gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes 
% Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
% Inputs:
% Fs        (sampling frequency associated with the data)-required
% nfft      (number of points in fft)-required
% fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
% Outputs:
% f         (frequencies)
% findx     (index of the frequencies in the full frequency grid). e.g.: If
% Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
% between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
% contain the indices in the frequency grid corresponding to frequencies <
% 100 Hz. In the case fpass=[0 500], findx=[1 512].
if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);
end

function [tapers,eigs]=dpsschk(tapers,N,Fs)
% Helper function to calculate tapers and, if precalculated tapers are supplied, 
% to check that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
% Usage: tapers=dpsschk(tapers,N,Fs)
% Inputs:
% tapers        (tapers in the form of: 
%                                   (i) precalculated tapers or,
%                                   (ii) [NW K] - time-bandwidth product, number of tapers) 
%
% N             (number of samples)
% Fs            (sampling frequency - this is required for nomalization of
%                                     tapers: we need tapers to be such
%                                     that integral of the square of each taper equals 1
%                                     dpss computes tapers such that the
%                                     SUM of squares equals 1 - so we need
%                                     to multiply the dpss computed tapers
%                                     by sqrt(Fs) to get the right
%                                     normalization)
% Outputs: 
% tapers        (calculated or precalculated tapers)
% eigs          (eigenvalues) 
if nargin < 3; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=dpss(N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
end

function J=mtfftc(data,tapers,nfft,Fs)
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if nargin < 4; error('Need all input arguments'); end;
data=change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
J=fft(data_proj,nfft)/Fs;   % fft of projected data
end

function [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,numsp1,numsp2)
% Function to compute lower and upper confidence intervals on the coherency 
% given the tapered fourier transforms, errchk, trialave.
%
% Usage: [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,numsp1,numsp2)
% Inputs:
% C     - coherence
% J1,J2 - tapered fourier transforms 
% err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
%                   p - p value for error estimates)
% trialave - 0: no averaging over trials/channels
%            1 : perform trial averaging
% numsp1    - number of spikes for data1. supply only if finite size corrections are required
% numsp2    - number of spikes for data2. supply only if finite size corrections are required
%
% Outputs: 
%          confC - confidence level for C - only for err(1)>=1
%          phistd - theoretical or jackknife standard deviation for phi for err(1)=1 and err(1)=2 
%                   respectively. returns zero if coherence is 1
%          Cerr  - Jacknife error bars for C  - only for err(1)=2
% Jackknife uses the following transform of the coherence
% z=sqrt(2*dim-2)atanh(C). Asymptotically (and for Gaussian data) var(z)=1.

if nargin < 5; error('Need at least 5 input arguments'); end;
if err(1)==0; error('Need err=[1 p] or [2 p] for error bar calculation'); end;
if nargout==4  && err(1)==1; error('Cerr contains Jackknife errors: only computed when err(1) is 2'); end;
[nf,K,Ch]=size(J1);
errchk=err(1);
p=err(2);
pp=1-p/2;
%
% Find the number of degrees of freedom
%
if trialave;
   dim=K*Ch;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   Ch=1;
   if nargin>=6 && ~isempty(numsp1) 
      totspikes1=sum(numsp1);
      dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
   end
   if nargin==7 && ~isempty(numsp2); 
      totspikes2=sum(numsp2);
      dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
   end;
   dof=min(dof1,dof2);
   J1=reshape(J1,nf,dim);
   J2=reshape(J2,nf,dim);
else
   dim=K;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   for ch=1:Ch;
      if nargin>=6 && ~isempty(numsp1);
         totspikes1=numsp1(ch); 
        dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
      end;
      if nargin==7 && ~isempty(numsp2);
         totspikes2=numsp2(ch);
        dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
      end;
      dof(ch)=min(dof1,dof2);
   end;
end;
%
% variance of the phase
%
%
% Old code is the next few lines - new code is in the if statement below
% beginning line 87
%
% if isempty(find((C-1).^2 < 10^-16));
%    phierr = sqrt((2./dof(ones(nf,1),:)).*(1./(C.^2) - 1));  
% else
%    phierr = zeros(nf,Ch);
% end  

%
% theoretical, asymptotic confidence level
%
if dof <= 2
   confC = 1;
else     
   df = 1./((dof/2)-1);
   confC = sqrt(1 - p.^df);
end;
%
% Phase standard deviation (theoretical and jackknife) and jackknife
% confidence intervals for C
%
if errchk==1;
   totnum=nf*Ch;
   phistd=zeros(totnum,1); 
   CC=reshape(C,[totnum,1]); 
   indx=find(abs(CC-1)>=1.e-16);
   dof=repmat(dof,[nf,1]);
   dof=reshape(dof,[totnum 1]); 
   phistd(indx)= sqrt((2./dof(indx).*(1./(C(indx).^2) - 1))); 
   phistd=reshape(phistd,[nf Ch]);
elseif errchk==2;
    tcrit=tinv(pp,dof-1);
    for k=1:dim; % dim is the number of 'independent' estimates
        indxk=setdiff(1:dim,k);
        J1k=J1(:,indxk,:);
        J2k=J2(:,indxk,:);
        eJ1k=squeeze(sum(J1k.*conj(J1k),2));
        eJ2k=squeeze(sum(J2k.*conj(J2k),2));
        eJ12k=squeeze(sum(conj(J1k).*J2k,2)); 
        Cxyk=eJ12k./sqrt(eJ1k.*eJ2k);
        absCxyk=abs(Cxyk);
        atanhCxyk(k,:,:)=sqrt(2*dim-2)*atanh(absCxyk); % 1-drop estimate of z
        phasefactorxyk(k,:,:)=Cxyk./absCxyk;
%         indxk=setdiff(1:dim,k);
%         J1jk=J1(:,indxk,:);
%         J2jk=J2(:,indxk,:);
%         eJ1jk=squeeze(sum(J1jk.*conj(J1jk),2));
%         eJ2jk=squeeze(sum(J2jk.*conj(J2jk),2));
%         eJ12jk=squeeze(sum(conj(J1jk).*J2jk,2)); 
%         atanhCxyjk(k,:,:)=sqrt(2*dim-2)*atanh(abs(eJ12jk)./sqrt(eJ1jk.*eJ2jk));
    end; 
    atanhC=sqrt(2*dim-2)*atanh(C); % z
    sigma12=sqrt(dim-1)*squeeze(std(atanhCxyk,1,1)); % Jackknife estimate std(z)=sqrt(dim-1)*std of 1-drop estimates
%     sigma12=sqrt(dim-1)*squeeze(std(atanhCxyjk,1,1));
    if Ch==1; sigma12=sigma12'; end;
    Cu=atanhC+tcrit(ones(nf,1),:).*sigma12;
    Cl=atanhC-tcrit(ones(nf,1),:).*sigma12;
    %Cerr(1,:,:) = tanh(Cl/sqrt(2*dim-2));
    Cerr(1,:,:) = max(tanh(Cl/sqrt(2*dim-2)),0); % This ensures that the lower confidence band remains positive
    Cerr(2,:,:) = tanh(Cu/sqrt(2*dim-2));
    %phistd=(2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk))));
    phistd = sqrt( (2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk)))) );
    if trialave; phistd=phistd'; end;
end
% ncrit=norminv(pp);
% phierr=zeros([2 size(phistd)]); 
% phierr(1,:,:)=phi-ncrit*phistd; phierr(2,:,:)=phi+ncrit*phistd;
end

end