import json
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
Parse the json file looking for a given process and build the list of files to process
"""
def getFilesForProcess(jsonUrl, tag, inDir):
    toReturn=[]
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()
    for proc in procList :
        for desc in proc[1] :
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                if dtag!=tag : continue
                split=getByLabel(d,'split',1)
                for segment in range(0,split) :
                    eventsFile=dtag
                    if split>1:
                        eventsFile=dtag + '_' + str(segment)
                    eventsFileUrl=inDir+'/'+eventsFile+'.root'
                    if(eventsFileUrl.find('/store/')==0)  :
                        eventsFileUrl = commands.getstatusoutput('cmsPfn ' + eventsFileUrl)[1]
                    toReturn.append(eventsFileUrl)
    return toReturn

"""
Parse the json file and retrieve all processes
"""
def getProcesses(jsonUrl,getMC=True,getData=False):
    toReturn=[]
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()
    for proc in procList :
        for desc in proc[1] :
            isData=getByLabel(desc,'isdata',False)
            if isData and not getData : continue
            if not isData and not getMC : continue
            data = desc['data']
            for d in data :
                toReturn.append( getByLabel(d,'dtag') )
    return toReturn
