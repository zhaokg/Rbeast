from . import Rbeast as cb
from numpy import ndarray, squeeze

def beast(Y, \
          start=1,
          deltat=1,
          season='harmonic',       # 'harmonic','dummy','svd (not supported yet)','none'
          freq  = float('nan'),
          scp_minmax = [0, 10] ,
          sorder_minmax = [0, 5],
          sseg_minlength = None, # an integer
          tcp_minmax = [0,10 ],
          torder_minmax = [0, 5],
          tseg_minlength = None,  # an integer
          detrend = False,
          deseasonalize = False,
          mcmc_seed  =0,
          mcmc_burbin = 200,
          mcmc_chains = 3,
          mcmc_thin = 5,
          mcmc_samples =8000,
          ci = False,
          precValue =1.5,
          precPriorType ='componentwise', #componentwise','uniform','constant','orderwise'
          print_options =True,
          print_progress = True,
          hasOutlier = False,
          ocp_max  = 10,
          gui     = False
        ):
      
      isNumpyInput = False;
      if   isinstance(Y, list) or isinstance(Y, tuple):
            isNumpyInput = False
      elif isinstance(Y, ndarray):
            isNumpyInput = True
      elif hasattr(Y,'to_numpy'):
            Y = getattr(Y,'to_numpy')()
            isNumpyInput = True
      else:
            raise ValueError('Unknown formats for the input Y.')
      
      if isNumpyInput:
            Y=squeeze(Y)
                 
    #......Start of displaying 'MetaData' ......
      metadata = lambda: None   ###Just get an empty object###
      metadata.isRegularOrdered = True
      metadata.season           = season
      metadata.startTime        = start
      metadata.deltaTime        = deltat
      #metadata.whichDimIsTime   = 1
      if (season != 'none'):
            metadata.period = deltat * freq
      metadata.missingValue     = float('nan')
      metadata.maxMissingRate   = 0.7500
      metadata.deseasonalize    = deseasonalize
      metadata.detrend          = detrend
      metadata.hasOutlierCmpnt  = hasOutlier
    #........End of displaying MetaData ........
    #......Start of displaying 'prior' ......
      prior = lambda: None   ###Just get an empty object###
      prior.modelPriorType	  = 1
      if season !='none' or season == None:
            prior.seasonMinOrder   = sorder_minmax[0]
            prior.seasonMaxOrder   = sorder_minmax[1]
            prior.seasonMinKnotNum = scp_minmax[0]
            prior.seasonMaxKnotNum = scp_minmax[1]
            prior.seasonMinSepDist = sseg_minlength
      prior.trendMinOrder	  = torder_minmax[0]
      prior.trendMaxOrder	  = torder_minmax[1]
      prior.trendMinKnotNum  = tcp_minmax[0]
      prior.trendMaxKnotNum  = tcp_minmax[1]
      prior.trendMinSepDist  = tseg_minlength
      prior.K_MAX            = 500
      prior.precValue        = precValue
      prior.precPriorType    = precPriorType
    #......End of displaying pripr ......
    #......Start of displaying 'mcmc' ......
      mcmc = lambda: None   ###Just get an empty object###
      mcmc.seed                      =  mcmc_seed
      mcmc.samples                   = mcmc_samples
      mcmc.thinningFactor            = mcmc_thin
      mcmc.burnin                    = mcmc_burbin
      mcmc.chainNumber               = mcmc_chains
      mcmc.maxMoveStepSize           = 6
      mcmc.trendResamplingOrderProb  = 0.1000
      mcmc.seasonResamplingOrderProb = 0.1700
      mcmc.credIntervalAlphaLevel    = 0.950
    #......End of displaying mcmc ......
    #......Start of displaying 'extra' ......
      extra = lambda: None   ###Just get an empty object###
      extra.dumpInputData        = True
      extra.whichOutputDimIsTime = 1
      extra.computeCredible      = True
      extra.fastCIComputation    = True
      extra.computeSeasonOrder   = True
      extra.computeTrendOrder    = True
      extra.computeSeasonChngpt  = True
      extra.computeTrendChngpt   = True
      extra.computeSeasonAmp     = False
      extra.computeTrendSlope    = True
      extra.tallyPosNegSeasonJump= False
      extra.tallyPosNegTrendJump = False
      extra.tallyIncDecTrendJump = False
      extra.printProgressBar     = True
      extra.printOptions         = True
      extra.consoleWidth         = 85
      extra.numThreadsPerCPU     = 2
      extra.numParThreads        = 0
    #......End of displaying extra ......
      if gui:
        o=cb.Rbeast('beastv4demo',Y, metadata, prior, mcmc, extra)
      else:
        o=cb.Rbeast('beastv4',Y, metadata, prior, mcmc, extra)      
      return (o)


