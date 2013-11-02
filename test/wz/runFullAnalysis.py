#!/usr/bin/env python
import os,sys
import json
import optparse
import commands

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal

"""
Loop over the inputs and launch jobs
"""
def runFullAnalysis(myExec, inDir, jsonUrl, params, outdir, onlyTag, queue) :

    if myExec.find('.py')>0 : myExec='python %s'%myExec
    
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()

    for proc in procList :
    
        for desc in proc[1] :

            isData=getByLabel(desc,'isdata',False)
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                if onlyTag!='' and dtag.find(onlyTag)<0 : continue
                split=getByLabel(d,'split',1)
                xsec = -1
                if not isData : xsec=getByLabel(d,'xsec',-1)

                for segment in range(0,split) :
                    if(split==1): 
                        eventsFile=inDir + '/' + dtag + '.root'
                    else:
                        eventsFile=inDir + '/' + dtag + '_' + str(segment) + '.root'

                    if(eventsFile.find('/store/')==0)  :
                        eventsFile = commands.getstatusoutput('cmsPfn ' + eventsFile)[1]

                    cmd='%s -i %s -o %s -x %f %s'%(myExec,eventsFile,outdir,xsec,params)
                    os.system(cmd)

def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-e', '--exe'        ,    dest='exe'                , help='Excecutable'                            , default=None)
    parser.add_option('-i', '--in'         ,    dest='inDir'              , help='Input directory'                        , default=None)
    parser.add_option('-s', '--sub'        ,    dest='queue'              , help='Batch queue (optional)'                 , default='')
    parser.add_option('-j', '--json'       ,    dest='json'               , help='A json file with the samples to analyze', default=None)
    parser.add_option('-o', '--out'        ,    dest='outdir'             , help='Output directory'                       , default='./results')
    parser.add_option('-t', '--tag'        ,    dest='onlyTag'            , help='Only matching'                          , default='')
    parser.add_option('-p', '--pars'       ,    dest='params'             , help='Extra parameters for the job'           , default='')
    (opt, args) = parser.parse_args()

    if opt.inDir is None or opt.exe is None or opt.json is None:
        parser.print_help()
        sys.exit(1)

    os.system('mkdir -p %s'%opt.outdir)
    runFullAnalysis(myExec=opt.exe, inDir=opt.inDir, jsonUrl=opt.json, params=opt.params, outdir=opt.outdir, onlyTag=opt.onlyTag, queue=opt.queue)

if __name__ == "__main__":
    main()
