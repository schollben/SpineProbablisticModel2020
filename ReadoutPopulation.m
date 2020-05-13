%ReadoutPopulation constructs a readout population object:
%POPOBJ = ReadoutPopulation(INPOP);
%
%   This object creates an optimal readout for an input population 
%   generated from InputPopulation. Assumptions of certain statistics 
%   about the input population: Gaussian with covariance that is invariant 
%   to the stimulus (does not change with stimulus).
%
%   The ReadoutPopulation object follows handle semantics; that is, methods
%   called on it affect the original object, not a copy of it. Also note
%   that ReadoutPopulation method names begin with a lowercase letter
%   (e.g., simulateResponses) while ReadoutPopulation property names begin
%   with an uppercase letter (e.g., Wts).
%
%   InputPopulation methods:
%       getStimId       -   Convert a real valued stimulus value into an
%                           index into the input population tuning curves.
%                           stimuli are treated as discrete categories in
%                           this population readout.
%       getExWtProb     -   Converts positive weights into a probability
%                           distribution
%       getInhWtProb    -   Converts negative weights into a probability
%                           distribution
%       getSynInPop     -   Uses excitatory weight probability to randomly
%                           generate a spine input population
%       plotWeights     -   Plot the optimal weights
%       getTuningCorr   -   Compute Pearson correlation between tuning of
%                           synaptic input population and a readout neuron  
%       getDeltaTuning  -   Get different in preference between input
%                           poulation and readout neuron
%       plotSynapticInputPopulation - plot tuning curves of input
%                                     population for a readout neuron
%       plotSomaTuning  -   Plot readout neuron stimulus tuning 
%       plotDeltaTuning -   Plot distribution of input preferences relative
%                           to readout neuron preference
%       doDecoding      -   Calculates readout accuracy, reported as 
%                           mean-square-error from INPOP via maximum a 
%                           posteriori estimation  
%
%   InputPopulation public fields:
%       Type            -   Type of readout (only Guassian support)
%                           Default: 'probgaussev'
%                           Force smooth weights: 'smoothgaussev' based on  
%                           Park and Pillow, PLoS Comp Bio (2011) 
%                           Input parameters ([p, d])
%       Wts             -   Readout weights for each neuron in the readout
%                           population
%       Offset          -   Baseline firing rate (DC term) for each neuron
%                           in the readout population
%       TuningCurves    -   The tuning curves evaluated at all possible
%                           stimuli
%
%
%   Example:
%      POPOBJ = ReadoutPopulation(INPOP);
%
%
%   Written by Jacob Yates 2018
%   Updated by Jacob Yates and Benjamin Scholl 2020

classdef ReadoutPopulation < handle
    
    properties
        Type
        Wts
        Offset
        TuningCurves
    end
    
    properties(Access=private)
        Ipop@InputPopulation
    end
    
    methods
        
        function obj = ReadoutPopulation(I, varargin)
            
            assert(isa(I, 'InputPopulation'), 'first argument must be an InputPopulation object')
            
            obj.Ipop = I;
            
            ip = inputParser();
            ip.KeepUnmatched = 1;
            ip.addParameter('Type', 'probgaussev', @(x) any(strcmp(x, {'probgaussev', 'smoothgaussev'})));
            ip.addParameter('asdParams', [0 0])
            ip.parse(varargin{:})
            
            obj.Type = ip.Results.Type;
            
            switch lower(obj.Type)
                case 'probgaussev'
                    % get optimal weights under assumption of normal
                    % distribution and equal variance
                    
                    obj.Wts = I.CovMat\I.TuningCurves';
                    obj.Offset = zeros(I.Params.NumStim,1);
                    
                    for iStim = 1:I.Params.NumStim
                        obj.Offset(iStim) = (-.5) * I.TuningCurves(iStim,:)*(I.CovMat\I.TuningCurves(iStim,:)') - log(I.Params.NumStim);
                    end
                    
                case 'smoothgaussev'
                    % here we can force smoothness in orientation space
                    % C_ij = exp(-p - (delta_ij / 2(d^2)))
                    % from Park and Pillow, PLoS Comp Bio (2011)
                    
                    params = ip.Results.asdParams;
                    
                    delta = abs(angle(exp(1i*bsxfun(@minus, I.mu, I.mu(:)))))*180/pi; % circular distance between two tuning curves
                    
                    CpriorInv = exp(-params(1)-(delta./(2*params(2)^2)));
                    CpriorInv = pinv(CpriorInv);
                    obj.Wts = (I.CovMat + CpriorInv)\I.TuningCurves';
                    obj.Offset = zeros(I.Params.NumStim,1);
                    
                    for iStim = 1:I.Params.NumStim
                        obj.Offset(iStim) = -.5 * I.TuningCurves(iStim,:)*(I.CovMat\I.TuningCurves(iStim,:)') - log(I.Params.NumStim);
                    end
                    
            end
            
            NumStim = obj.Ipop.Params.NumStim;
            % set 1 readout neuron per stimulus category %
            NumNeurons = obj.Ipop.Params.NumStim;
            obj.TuningCurves = zeros(NumStim, NumNeurons);
            
            for j = 1:NumStim
                
                ak = obj.Ipop.TuningCurves(j,:)*obj.Wts + obj.Offset';
                
                for i = 1:NumNeurons
                    
                    pp = exp(ak(i))/sum(exp(ak));
                    
                    obj.TuningCurves(j,i) = pp;
                    
                end
            end
            
            
        end % constructor
       
        function id = getStimId(obj, stim)
            [~, id] = min(abs(obj.Ipop.circDist(stim, obj.Ipop.Stim(:))));
        end % getStimId
        
        
        function wprob = getExWtProb(obj, stim)
            if nargin < 2
                stim = 0;
            end
            
            id = obj.getStimId(stim);
            
            wprob = obj.Wts(:,id);
            wprob(wprob < 0) = 0;
            wprob = wprob/sum(wprob);
            
        end % getExcitatoryWeighProbabilities
        
        
        function wprob = getInhWtProb(obj, stim)
            if nargin < 2
                stim = 0;
            end
            
            id = obj.getStimId(stim);
            
            wprob = obj.Wts(:,id);
            wprob(wprob > 0) = 0;
            wprob = wprob/sum(wprob);
            
        end % getInhibitoryWeighProbabilities
        
        
        function [sip, inputId] = getSynInPop(obj, nMeasuredSpines, stim)
            if nargin < 3
                stim = 0;
                if nargin < 2
                    nMeasuredSpines = 100;
                end
            end
            
            wprob = obj.getExWtProb(stim);
            
            inputId = randsample(1:obj.Ipop.NumNeurons, nMeasuredSpines, true, wprob);
            
            sip = obj.Ipop.TuningCurves(:, inputId);
            
            sip = sip + randn(size(sip))*.1*max(obj.Ipop.TuningCurves(:));
            
        end % getSynInPop
        
        
        function [rho, sip, tc, inputId] = getTuningCorr(obj, nMeasuredSpines, stim, subsample)
            % [rho, sip, tc, inputId] = getTuningCorr(obj, nMeasuredSpines, stim, subsample)
            % get the tuning correlation between the synaptic input
            % population and the soma of a neuron tuned for stim
            if nargin < 4
                subsample = 1;
                if nargin < 3
                    stim = 0;
                    if nargin < 2
                        nMeasuredSpines = 100;
                    end
                end
            end
            
            [sip, inputId] = obj.getSynInPop(nMeasuredSpines, stim);
            sip = bsxfun(@rdivide, sip, sum(sip));
            
            % get the neuron tuned for this stimulus
            stimId = obj.getStimId(stim);
            tc = obj.TuningCurves(:,stimId);
            tc = tc/sum(tc);
            
            sip = sip(1:subsample:end,:);
            tc = tc(1:subsample:end);
            rho = corr(sip, tc,'type','Pearson');
            
        end
        
        
        function ds = getDeltaTuning(obj, nMeasuredSpines, stim)
            
            % get the tuning correlation between the synaptic input
            % population and the soma of a neuron tuned for stim
            if nargin < 3
                stim = 0;
                if nargin < 2
                    nMeasuredSpines = 100;
                end
            end
            
            [~, inputId] = obj.getSynInPop(nMeasuredSpines, stim);
            
            ds = obj.Ipop.circDist(obj.Ipop.mu(inputId), stim);
            
        end % getDeltaTuning
        
        
        function plotHandles = plotWeights(obj, stim)
            
            if ~exist('stim', 'var')
                stimId = 1:numel(obj.Ipop.Stim);
            else
                stimId = obj.getStimId(stim);
            end
            
            plotHandles = plot(obj.Ipop.mu, obj.Wts(:,stimId)');
            
        end % plotWeights
        
        
        function plotWeightsSelectivity(obj,I)
            
            tc = I.TuningCurves;
            angs = I.Stim + pi/2;
            vecStrength = [];
            for j = 1:size(tc,2)
                amp = tc(:,j)';
                vecStrength(j) = sqrt(sum(cos(angs).*amp)^2 + sum(sin(angs).*amp)^2)/sum(amp);
            end%vector strength
            
            nStim = length(I.Stim);
            numNeurons = I.NumNeurons;
            
            rr = zeros(nStim*numNeurons,1);
            ds = zeros(nStim*numNeurons,1);
            wprob = zeros(nStim*numNeurons,1);
            
            for n = 1:nStim
                inds = (1:numNeurons) + numNeurons*(n - 1);
                rr(inds) = vecStrength;
                ds(inds) = abs(I.circDist(I.mu, I.Stim(n)));          
                wprob(inds) = obj.getExWtProb(I.Stim(n));
            end
            ds(ds>pi/2) = ds(ds>pi/2) - pi/2;
            
            figure(101); clf
            subplot(1,2,1)
            plot(rr(wprob>0),wprob(wprob>0),'.k')
            subplot(1,2,2)
            plot(ds(wprob>0),wprob(wprob>0),'.k')
            shg
        end % plotWeights against Input Population seletivity/preference
        
        
        function plotSynapticInputPopulation(obj, NumMeasuredSpines, stim, version)
            if nargin < 4
                version = 1;
                if nargin < 3
                    stim = 0;
                    if nargin < 2
                        NumMeasuredSpines = 100;
                    end
                end
            end
            
            sip = obj.getSynInPop(NumMeasuredSpines, stim);
            switch version
                case 1
                    plot(obj.Ipop.Stim, sip)
                    axis tight
                    xlabel('\theta')
                    ylabel('Response')
                case 2
                    [~, mid] = max(sip);
                    [~, id] = sort(mid);
                    imagesc(obj.Ipop.Stim,[], sip(:,id)')
                    xlabel('\theta')
            end
            title('Synaptic Input Population', 'Fontweight', 'Normal')
            
        end % plotSynapticInputPopulation
        
        
        function plotHandles = plotSomaTuning(obj, stim)
            
            if nargin < 2
                stim = 0;
            end
            stimId = obj.getStimId(stim);
            tc = obj.TuningCurves(:,stimId);
            plotHandles = plot(obj.Ipop.Stim, tc);
            xlabel('\theta')
            ylabel('Response')
            title('Soma Tuning Curve', 'Fontweight', 'Normal')
        end % plotSomaTuning
        
        
        function [cnt, bins] = plotDeltaTuning(obj, nMeasuredSpines, stim, binSize)
            
            if nargin < 4
                binSize = pi/20;
                if nargin < 3
                    stim = 0;
                    if nargin < 2
                        nMeasuredSpines = 100;
                    end
                end
            end
            
            ds = obj.getDeltaTuning(nMeasuredSpines, stim);
            
            binEdges = 0:binSize:max(ds)+binSize;
            cnt = histc(abs(ds), binEdges);
            bins = binEdges + binSize/2;
            
            stairs(bins, cnt);
            
        end % plotDeltaTuning
        
        
        function [MSE] = doDecoding(obj,I, varargin)
            % Decode stimulus value from the response of input population
           
            numRepeats = 100;
            
            R  = [];
            id = [];
            for stimId = 1:numel(I.Stim)
                id = [id; stimId*ones(numRepeats,1)];
                R = [R; I.simulateResponses(stimId, numRepeats)];
            end
            wproj = R*obj.Wts + obj.Offset';
            [~, ihat] = max(wproj, [],2); %MAP estimation
            mu = I.Stim(ihat); 
            MSE = nanmean(angle(exp(1i*(mu(:)-I.Stim(id)'))).^2);%accuracy
            
        end % doDecoding
        
        
    end % methods
    
end % classdef