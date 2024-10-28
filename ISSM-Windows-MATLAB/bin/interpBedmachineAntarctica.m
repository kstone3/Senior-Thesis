function output = interpBedmachineAntarctica(X,Y,string,method,ncdate)
%INTERPBEDMACHINEANTARCTICA - interpolate BedMachine data onto X and Y
%
%   Examples:
%      bed       = interpBedmachineAntarctica(X,Y,'bed');
%      surface   = interpBedmachineAntarctica(X,Y,'surface');
%      thickness = interpBedmachineAntarctica(X,Y,'thickness');
%      mask      = interpBedmachineAntarctica(X,Y,'mask');
%      mask      = interpBedmachineAntarctica(X,Y,'mask','nearest','../Data/BedMachineAntarctica_2020-07-15_v02.nc');
%
%   - mask:   0 ocean, 1 land (ice free), 2 grounded ice, 3 floating ice
%   - source: 1 IBCSO/RTopo-2, 2 MC, 3 interpolation, 4 hydrostatic eq, 
%             5 Streamline diffusion, 6 Gravity inversion
%   - optional 4th input argument: interpolation method.
%             Supported interpolation methos: 'linear','cubic','nearest'
%   - optional 5th input argument: path to dataset.
%
% Version 11/30/2018 Mathieu Morlighem mmorligh@uci.edu

if nargin<3, string = 'bed'; end
if nargin<4
	if strcmp(string,'mask') | strcmp(string,'source')
		method='nearest'; % default method
	else
		method='cubic'; % default method
	end
end
if nargin<5
	ncdate='2020-07-15'; %BedMachine v2
	ncdate='v3.5';       %Official v3 release
end
basename = 'BedMachineAntarctica';

if nargin==5
	ncfile = ncdate;
else
	%List of common paths to try
	paths = {...
		['/u/astrid-r1b/ModelData/BedMachine/' basename '-' ncdate '.nc'],...
		['/home/ModelData/Antarctica/BedMachine/' basename '-' ncdate '.nc'],...
		['/totten_1/ModelData/Antarctica/BedMachine/' basename '-' ncdate '.nc'],...
		['/Users/larour/ModelData/BedMachine/' basename '-' ncdate '.nc'],...
		['./' basename '-' ncdate '.nc'],...
		};

	found = 0;
	for i=1:numel(paths)
		if exist(paths{i},'file')
			ncfile = paths{i};
			found = 1;
			break;
		end
	end

	if ~found
		error(['Could not find ' basename '-' ncdate '.nc, you can add the path to the list or provide its path as a 5th argument']);
	end
end

disp(['   -- BedMachine Antarctica version: ' ncdate]);
xdata = double(ncread(ncfile,'x'));
ydata = double(ncread(ncfile,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
if isempty(posx), posx=numel(xdata); end
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
if isempty(posy), posy=numel(ydata); end
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

if strcmp(string,'icemask'),
	disp(['   -- BedMachine Antarctica: loading ' string]);
	%data  = double(ncread(ncfile,'mask'))';
	data  = double(ncread(ncfile,'mask',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
	xdata=xdata(id1x:id2x);
	ydata=ydata(id1y:id2y);
	%ice ocean interface is between 0 and 3, so we might get some 1 by interpolating
	data(find(data==3))=0;
else
	disp(['   -- BedMachine Antarctica: loading ' string]);
	%data  = double(ncread(ncfile,string))';
	data  = double(ncread(ncfile,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
	xdata=xdata(id1x:id2x);
	ydata=ydata(id1y:id2y);
end

disp(['   -- BedMachine Antarctica: interpolating ' string]);
disp(['       -- Interpolation method: ' method]);
if strcmp(string,'mask') | strcmp(string,'source'),
	%Need nearest neighbor to avoid interpolation between 0 and 2
	output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');
	%tic
	%output = FastInterp(xdata,ydata,data,X,Y,'nearest');
	%toc
else
	%disp('InterpFromGrid');
	%tic
	%output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'cubic'); 
	output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),method); % now the interpolation method can be defined by the user
	%toc
	%disp('FastInterp');
	%tic
	%output = FastInterp(xdata,ydata,data,X,Y,'bilinear');
	%toc
end

end
function zi = FastInterp(x,y,data,xi,yi,method)

	%get data size
	[M N] = size(data);

	% Get X and Y library array spacing
	ndx = 1/(x(2)-x(1));    ndy = 1/(y(2)-y(1));
	% Begin mapping xi and yi vectors onto index space by subtracting library
	% array minima and scaling to index spacing

	xi = (xi - x(1))*ndx;       yi = (yi - y(1))*ndy;

	% Fill Zi with NaNs
	zi = NaN(size(xi));

	if strcmpi(method,'nearest'),
		% Find the nearest point in index space
		rxi = round(xi)+1;  ryi = round(yi)+1;
		% Find points that are in X,Y range
		flag = rxi>0 & rxi<=N & ~isnan(rxi) & ryi>0 & ryi<=M & ~isnan(ryi);
		% Map subscripts to indices
		ind = ryi + M*(rxi-1);
		zi(flag) = data(ind(flag));

	else %Bilinear

		% Transform to unit square
		fxi = floor(xi)+1;  fyi = floor(yi)+1; % x_i and y_i
		dfxi = xi-fxi+1;    dfyi = yi-fyi+1;   % Location in unit square

		% flagIn determines whether the requested location is inside of the data arrays
		flagIn = fxi>0 & fxi<N & ~isnan(fxi) & fyi>0 & fyi<M & ~isnan(fyi);

		%Toss all out-of-bounds variables now to save time
		fxi  = fxi(flagIn);  fyi  = fyi(flagIn);
		dfxi = dfxi(flagIn); dfyi = dfyi(flagIn);

		%Find bounding vertices
		ind1 = fyi + M*(fxi-1);     % indices of (  x_i  ,  y_i  )
		ind2 = fyi + M*fxi;         % indices of ( x_i+1 ,  y_i  )
		ind3 = fyi + 1 + M*fxi;     % indices of ( x_i+1 , y_i+1 )
		ind4 = fyi + 1 + M*(fxi-1); % indices of (  x_i  , y_i+1 )

		% Bilinear interpolation
		zi(flagIn) = ...
			data(ind1).*(1-dfxi).*(1-dfyi) + ...
			data(ind2).*dfxi.*(1-dfyi) + ...
			data(ind4).*(1-dfxi).*dfyi + ...
			data(ind3).*dfxi.*dfyi;
	end
end
