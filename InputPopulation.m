%InputPopulation constructs an input population object:
%INPOP = InputPopulation();
%
%   This object supports the simulation of a tuning curve model of a neural
%   population. The constructor specifies the tuning curve
%   parameterization, the correlation structure, homegeneity of the
%   population and size of the population.
%
%   The InputPopulation object follows handle semantics; that is, methods 
%   called on it affect the original object, not a copy of it. Also note that
%   InputPopulation method names begin with a lowercase letter (e.g., plotTuning)
%   while InputPopulation property names begin with an uppercase letter 
%   (e.g., NumNeurons).
%
%   InputPopulation methods:
%       plotTuning        - plot the Input population tuning curves
%       simulateResponses - simulate Poisson spiking activity across stimulus trials    
%
%   InputPopulation public fields:
%       NumNeurons      -   Number of neurons in the population
%       Homogenous      -   (Boolean) Specify homogenous or hetergenous 
%                           tunining curves
%       StimSupport     -   Stimulus range in radians
%       Params          -   Structure that contains information about the
%                           construction of this specific object.
%       TuningFun       -   Tuning curve function
%       TuningParam     -   Tuning curve parameters
%       CorrMat         -   The correlation matrix of the population used
%                           to construct the covariance
%       TuningCurves    -   The tuning curves evaluated at all possible
%                           stimuli
%       maxNoiseCorr    -   Maximum amplitude of limited-range pariwise
%                           correlations between neurons
%   Example:
%       INPOP = InputPopulation(...
%           'NumNeurons',       1000, ...
%           'Homogenous',       true, ...
%           'NumStim',          90, ...
%           'StimSupport',      [-pi/2, pi/2], ...
%           'Tuning',           'Von Mises', ...
%           'maxNoiseCorr',     0.2); 
%
%       
%   Written by Jacob Yates 2018
%   Updated by Jacob Yates and Benjamin Scholl 2020

classdef InputPopulation < handle
    
    properties
        NumNeurons
        Params
        CorrMat
        CovMat
        TuningCurves
        Stim
        mu
        sigma
        FixedParamsValues
        InputPopVals
    end
    
    properties(Access=private)
       alpha
       beta
       TuningFun
       TuningParams
    end
    
    
    methods
        % --- constructor: build the input population
        function obj = InputPopulation(varargin)
            
            ip = inputParser();
            ip.KeepUnmatched = 1;
            ip.addParameter('NumNeurons', 100, @isnumeric);
            
            % What type of tuning curves? These are 1D tuning curves currently
            acceptedTuningTypes = {'Gaussian', 'Von Mises', 'Orientation'};
            ip.addParameter('Tuning', 'Von Mises', @(x) any(strcmpi(x, acceptedTuningTypes)));
            
            % What correlation structure does the population have?
            acceptedCorrelationStructures = {'none', 'limited-range'};
            ip.addParameter('Correlations', 'limited-range', @(x) any(strcmpi(x, acceptedCorrelationStructures)));
            ip.addParameter('Homogenous', false)
            ip.addParameter('FixTuningParams', false)
            ip.addParameter('FixedParamsValues', [5 0.3], @isnumeric);

            % what is the noise type?
            acceptedNoiseTypes = {'Poisson-like'};
            ip.addParameter('Noise', 'Poisson-like', @(x) any(strcmpi(x, acceptedNoiseTypes)));
            
            % setup the stimulus support
            ip.addParameter('NumStim', 32)
            ip.addParameter('StimSupport', [-pi/2 pi/2])
            
            ip.parse(varargin{:})
            
            obj.NumNeurons = ip.Results.NumNeurons;
            obj.Params = ip.Results;
            obj.FixedParamsValues = ip.Results.FixedParamsValues;
            
            switch lower(obj.Params.Tuning)
                                
                case 'gaussian'
                    obj.TuningFun = @(s, mu, alpha, beta, sigma) alpha + beta .* exp ( - (s(:) - mu).^2 ./ 2 ./sigma.^2 );
                    
                case 'von mises'
                    obj.TuningFun = @(s, mu, alpha, beta, sigma) bsxfun(@plus, alpha , bsxfun(@times, beta , exp ( bsxfun(@times, sigma,(cos(bsxfun(@minus, s(:),mu)) - 1)))));
                   
                case 'orientation'
                    obj.TuningFun = @(s, mu, alpha, beta, sigma) bsxfun(@plus, alpha, bsxfun(@times, beta , exp ( bsxfun(@times, sigma, (cos(bsxfun(@minus, s(:), mu)).^2 - 1)))));
            end
            
            % -- setup tuning curve parameters
            fields = fieldnames(ip.Unmatched);
            tuningArgs = cell(1,numel(fields)*2);
            for iField = 1:numel(fields)
                tuningArgs{(iField - 1) * 2 + 1} = fields{iField};
                tuningArgs{(iField - 1) * 2 + 2} = ip.Unmatched.(fields{iField});
            end
            
            ip2 = inputParser();
            ip2.KeepUnmatched = 1;
            ip2.addParameter('minFiringRate', 0, @(x) x >=0);
            ip2.addParameter('firingRateScale', 5, @(x) x >0);
            ip2.addParameter('bandwidthScale', 2, @(x) x>0);
            ip2.addParameter('bandwidthMin', 1, @(x) x>=0);
            ip2.parse(tuningArgs{:});
            
            obj.TuningParams = ip2.Results;
            
            if obj.Params.Homogenous
                % homogenous input population case
                obj.alpha = obj.TuningParams.minFiringRate + ones(1, obj.NumNeurons);
                obj.beta  = obj.TuningParams.firingRateScale + ones(1, obj.NumNeurons) + 1;
                obj.sigma = obj.TuningParams.bandwidthScale * ones(1, obj.NumNeurons) + obj.TuningParams.bandwidthMin;

                obj.mu = linspace(obj.Params.StimSupport(1), obj.Params.StimSupport(2), obj.Params.NumNeurons);
            else
                % heterogenous input population case
                obj.alpha = obj.TuningParams.minFiringRate + rand(1, obj.NumNeurons);
                obj.beta  = obj.TuningParams.firingRateScale * rand(1, obj.NumNeurons) + 1;
                
                % tuning curve widths drawn from a lognormal distribution
                % matched to Ringach et al., 2002, [-1 0.3]
                % sample Half Bandwidth at 1/sqrt(2) height (Fig. 2)
                mu_ = -1;
                sig_ = 0.3;     
                
                if obj.Params.FixTuningParams
                    % using defined paramteres here
                    tuningWidthVar = obj.FixedParamsValues(2);
                    sig_ = tuningWidthVar;
                    
                    firingRateVar = obj.FixedParamsValues(1);
                    C = (obj.TuningParams.firingRateScale-firingRateVar)*.5;
                    C(C<0) = 0;
                    obj.beta  = firingRateVar*rand(1, obj.NumNeurons) + 1 + C;
                end
                
                HBW = lognrnd(mu_, sig_, [1, obj.NumNeurons]); % units in radians
                %%%% in case you want to check that this distribution matches
                %%%% plot( (0:90), lognpdf((0:90)/90*(pi/2), -1 , .3))
                
                % keep standard deviation of tuning population amplitudes and widths (in deg)
                obj.InputPopVals = [range(obj.beta) sig_*(180/pi)];

                % convert half bandwidth (HBW) to the sigma parameter in our tuning curve model
                % Ringach et al. (2002) report bandwith at 1/sqrt(2) of max
                switch lower(obj.Params.Tuning)
                    
                    case 'gaussian'
                        hw2k = @(hw) hw./(2*sqrt(log(2)));
                        
                    case 'von mises'
                        hw2k = @(hw) -log(sqrt(2))./(cos(hw) - 1);
                        
                    case 'orientation'
                        hw2k = @(hw) -log(sqrt(2))./(cos(hw).^2 - 1);
                end
                
                obj.sigma = hw2k(HBW);
                
                obj.mu = linspace(obj.Params.StimSupport(1), obj.Params.StimSupport(2), obj.Params.NumNeurons);
                
            end
           
            
            % --- setup stimulus support
            obj.Stim = linspace(obj.Params.StimSupport(1), obj.Params.StimSupport(2), obj.Params.NumStim);
            
            obj.TuningCurves = obj.TuningFun(obj.Stim, obj.mu, obj.alpha, obj.beta, obj.sigma);
           

            % --- setup noise correlations
            ip3 = inputParser();
            ip3.KeepUnmatched = 1;
            ip3.addParameter('maxNoiseCorr', 0.2);
            ip3.addParameter('eps', .01);
            ip3.parse(tuningArgs{:});
            
            obj.Params.maxNoiseCorr = ip3.Results.maxNoiseCorr;
            
            switch lower(obj.Params.Correlations)
                case 'none'
                    % diagonal correlation matrix (no correlations)
                    obj.CorrMat = eye(obj.NumNeurons);
                case 'limited-range'
                    % correlation matrix
                    obj.CorrMat = obj.Params.maxNoiseCorr*exp(-abs(obj.circDist(obj.mu, obj.mu(:))));
                    obj.CorrMat = obj.CorrMat + (1 - obj.Params.maxNoiseCorr) * eye(obj.NumNeurons);
            end
            
            % build covariance matrix
            switch lower(obj.Params.Noise)
                case 'poisson-like'
                    poisson_std_dev = sqrt(diag(mean(obj.TuningCurves))); % std = sqrt(mean) for poisson
                    obj.CovMat = poisson_std_dev*obj.CorrMat*poisson_std_dev;
            end
            
            
        end % constructor
        
        function R = simulateResponses(obj, stimId, nTrials)
            
            if nargin < 3
                nTrials = 1;
            end
            
            mus = obj.TuningCurves(stimId,:);
            poisson_std_dev = sqrt(diag(mus)); % std = sqrt(mean) for poisson
            
            C = poisson_std_dev*obj.CorrMat*poisson_std_dev;
            
            R = mvnrnd(mus, C, nTrials);
        end
        
        function plotTuning(obj)
            
            plot(obj.Stim, obj.TuningCurves)
            axis tight
            xlabel('\theta')
            ylabel('Firing Rate')
            title('Population', 'FontWeight', 'normal')
            
        end % plotTuning
        
        function plotTuningSparse(obj)
            
            plot(obj.Stim, obj.TuningCurves(:,1:50:end),'color',[.5 .5 .5])
            axis tight
            xlabel('\theta')
            ylabel('Firing Rate')
            title('Population', 'FontWeight', 'normal')
            
        end % plotTuning
        
    end % methods
    
    
    methods(Static)
        
        function d = circDist(x,y)
            
            d = angle(exp(1i * bsxfun(@minus, x, y)));
            
        end
        
        function d = circDistDeg(x,y)
            
            x = x/180*pi;
            y = y/180*pi;
            
            d = circDist(x,y)/pi*180;
        end
        
    end
    
end % classdef

