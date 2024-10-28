%CFDRAGCOEFFABSGRAD class definition
%
%   Usage:
%      cfdragcoeffabsgrad=cfdragcoeffabsgrad();
%      cfdragcoeffabsgrad=cfdragcoeffabsgrad('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%							'model_string','Surface',...
%                    'weights',ones(md.mesh.numberofvertices,1),...
%                    'weights_string','WeightsSurfaceObservations');
%
%

classdef cfdragcoeffabsgrad
	properties (SetAccess=public)
		%cfdragcoeffabsgrad
		name               = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		weights            = NaN; %weight coefficients for every vertex
		weights_string     = ''; %string to identify this particular set of weights
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
		end % }}}
		function self = cfdragcoeffabsgrad(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.weights=getfieldvalue(options,'weights',NaN);
				self.weights_string=getfieldvalue(options,'weights_string','');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('cfdragcoeffabsgrad error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);

			md = checkfield(md,'fieldname','self.weights','field',self.weights,'timeseries',1,'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   TimeMisfit:\n'));

			fielddisplay(self,'name','identifier for this cfdragcoeffabsgrad response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'weights','weights (at vertices) to apply to the cfdragcoeffabsgrad');
			fielddisplay(self,'weights_string','string for weights for identification purposes');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

		WriteData(fid,prefix,'data',self.name,'name','md.cfdragcoeffabsgrad.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.cfdragcoeffabsgrad.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.weights,'name','md.cfdragcoeffabsgrad.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.weights_string,'name','md.cfdragcoeffabsgrad.weights_string','format','String');
		end % }}}
	end
end
