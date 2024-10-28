%MISFIT class definition
%
%   Usage:
%      cfsurfacesquare=cfsurfacesquare();
%      cfsurfacesquare=cfsurfacesquare('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%							'model_string','Surface',...
%                    'observation_string','SurfaceObservations',...
%                    'observation',md.geometry.surface,...
%                    'weights',ones(md.mesh.numberofvertices,1),...
%                    'weights_string','WeightsSurfaceObservations',...
%							'datatime',time);
%
%

classdef cfsurfacesquare
	properties (SetAccess=public)
		%cfsurfacesquare
		name               = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		model_string       = ''; %string for field that is modeled
		observation        = NaN; %observed field that we compare the model against
		observation_string = ''; %string for observed field.
		weights            = NaN; %weight coefficients for every vertex
		weights_string     = ''; %string to identify this particular set of weights
		datatime				 = 0; %time in years from start that the data is from 
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
			if ~isnan(self.observation)
				self.observation=project3d(md,'vector',self.observation,'type','node');
			end
		end % }}}
		function self = cfsurfacesquare(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.model_string=getfieldvalue(options,'model_string');
				self.observation=getfieldvalue(options,'observation',NaN);
				self.observation_string=getfieldvalue(options,'observation_string');
				self.weights=getfieldvalue(options,'weights',NaN);
				self.weights_string=getfieldvalue(options,'weights_string','');
				self.datatime = getfieldvalue(options, 'datatime');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.datatime = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('cfsurfacesquare error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);

			md = checkfield(md,'fieldname','self.observation','field',self.observation,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','self.weights','field',self.weights,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','self.datatime','field',self.datatime,'<=',md.timestepping.final_time);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   TimeMisfit:\n'));

			fielddisplay(self,'name','identifier for this cfsurfacesquare response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'model_string','string for field that is modeled');
			fielddisplay(self,'observation','observed field that we compare the model against');
			fielddisplay(self,'observation_string','observation string');
			fielddisplay(self,'weights','weights (at vertices) to apply to the cfsurfacesquare');
			fielddisplay(self,'weights_string','string for weights for identification purposes');
			fielddisplay(self,'datatime','time to compare data to model for misfit');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

		WriteData(fid,prefix,'data',self.name,'name','md.cfsurfacesquare.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.cfsurfacesquare.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.model_string,'name','md.cfsurfacesquare.model_string','format','String');
		WriteData(fid,prefix,'data',self.observation,'name','md.cfsurfacesquare.observation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.observation_string,'name','md.cfsurfacesquare.observation_string','format','String');
		WriteData(fid,prefix,'data',self.weights,'name','md.cfsurfacesquare.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.weights_string,'name','md.cfsurfacesquare.weights_string','format','String');
		WriteData(fid,prefix,'data',round(self.datatime*md.constants.yts),'name','md.cfsurfacesquare.datatime','format','Double');
		end % }}}
	end
end
