import obspy as obs
from obspy.geodetics.base import gps2dist_azimuth as gps
from obspy.core.event import Catalog
from tqdm import tqdm
from json import load
import warnings
warnings.filterwarnings("ignore")

with open("config.json") as f:
    configs = load(f)

catRef = obs.read_events("cat1.dat")
catCom = obs.read_events("cat2.dat")
catTar = Catalog()
catComEvents = Catalog()
timeShift = configs["TimeShift"] # acceptable time shift between two common events (s)
epicShift = configs["EpicentralShift"] #  acceptable epicenter shift between two common events (km)

def computeDiff(RefEvent, ComEvent):
    """Compute time and epicenter difference between "reference" and "compared" events,
       and decide whether two events are common or not.

    Args:
        RefEvent (obspy.event): reference event
        ComEvent (obspy.event): compared event

    Returns:
        bool: true if two events are common
    """
    RefEventPreferredOrigin = RefEvent.preferred_origin()
    ComEventPreferredOrigin = ComEvent.preferred_origin()    
    dT = abs(RefEventPreferredOrigin.time - ComEventPreferredOrigin.time)
    lat1 = RefEventPreferredOrigin.latitude
    lon1 = RefEventPreferredOrigin.longitude
    lat2 = ComEventPreferredOrigin.latitude
    lon2 = ComEventPreferredOrigin.longitude
    dD = abs(gps(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2)[0]*1e-3)
    if dT < timeShift and dD < epicShift:
        return True
    else:
        return False

def updateEvent(RefEvent, ComEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases):
    """Update "reference" event with new phases comes from "compared" event

    Args:
        RefEvent (obspy.event): "reference" event
        ComEvent (obspy.event): "compared" event
        numNewPhaseP (int): number of new P phases used to update "reference" event
        numNewPhaseS (int): number of new S phases used to update "reference" event
        numRegardedPhases (int): number of phases not used to update "reference" event

    Returns:
        obspy.event: updated "reference" event.
    """
    RefEventPreferredOrigin = RefEvent.preferred_origin()
    ComEventPreferredOrigin = ComEvent.preferred_origin()
    RefEventPicks = RefEvent.picks
    ComEventPicks = ComEvent.picks
    RefEventArrivals = RefEventPreferredOrigin.arrivals
    ComEventArrivals = ComEventPreferredOrigin.arrivals
    RefEventPickRemarks = ["_".join([pick.waveform_id.station_code, pick.phase_hint]) for pick in RefEventPicks]
    ComEventPickRemarks = ["_".join([pick.waveform_id.station_code, pick.phase_hint]) for pick in ComEventPicks]
    newPicks = list(set(ComEventPickRemarks) - set(RefEventPickRemarks))
    numRegardedPhases += len(ComEventPickRemarks) - len(newPicks)
    newPicksIndex = [ComEventPickRemarks.index(newPick) for newPick in newPicks]
    for newPickIndex in newPicksIndex:
        newPickResource_id = ComEventPicks[newPickIndex].resource_id
        newArrivals = [arrival for arrival in ComEventArrivals if arrival.pick_id == newPickResource_id]
        RefEventPicks.append(ComEventPicks[newPickIndex])
        RefEventArrivals.extend(newArrivals)
        if "P" in ComEventPicks[newPickIndex].phase_hint:
            numNewPhaseP += 1
        elif "S" in ComEventPicks[newPickIndex].phase_hint:
            numNewPhaseS += 1            
    RefEvent.picks = RefEventPicks
    RefEvent.preferred_origin().arrivals = RefEventArrivals
    return RefEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases

# starts application
numCommonEvents = 0
numNewPhaseP = 0
numNewPhaseS = 0
numRegardedPhases = 0

print("+++ Starts merging catalogs...")
for eventID, RefEvent in enumerate(tqdm(catRef)):
    for ComEvent in catCom:
        if computeDiff(RefEvent, ComEvent):
            RefEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases = updateEvent(RefEvent, ComEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases)
            catComEvents.append(RefEvent)
            numCommonEvents += 1
            break
    catTar.append(RefEvent)

catTar.write("updatedCatalog.dat", format="NORDIC")
if configs["OutputCommonEventsCatalog"]:
    catComEvents.write("commonEventsCatalog.dat", format="NORDIC")

# summary
print("+++ Number of common events: {numCommonEvents}".format(numCommonEvents=numCommonEvents))
print("+++ Number of new P phases: {numNewPhaseP}".format(numNewPhaseP=numNewPhaseP))
print("+++ Number of new S phases: {numNewPhaseS}".format(numNewPhaseS=numNewPhaseS))
print("+++ Number of regarded phases: {numRegardedPhases}".format(numRegardedPhases=numRegardedPhases))
