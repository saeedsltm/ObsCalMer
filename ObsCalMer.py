import warnings
from json import load

import obspy as obs
from obspy.core.event import Catalog
from obspy.geodetics.base import gps2dist_azimuth as gps
from tqdm import tqdm

warnings.filterwarnings("ignore")

# read config file
with open("config.json") as f:
    configs = load(f)

RefCat = obs.read_events(configs["ReferenceCatalogFileName"])
ComCat = obs.read_events(configs["ComparedCatalogFileName"])
TarCat = Catalog()
ComEventsCat = Catalog()
TimeShift = configs["TimeShift"]
EpicShift = configs["EpicentralShift"]


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
    if dT < TimeShift and dD < EpicShift:
        return True
    else:
        return False


def readExtra(pick):
    """reads pick weight 

    Args:
        pick (event.pick): an obspy event pick

    Returns:
        str: pick weight 
    """

    try:
        weight = pick.extra.get("nordic_pick_weight", "0")
        if isinstance(weight, type({})):
            weight = weight["value"]
    except AttributeError:
        weight = "0"
    return weight


def manageNewPicks(RefEventPicks, ComEventPicks, numRegardedPhases):
    """Manage new picks

    Args:
        RefEventPicks (event.picks): an obspy event picks, for "reference" event,
        ComEventPicks (event.picks): an obspy event picks, for "compared" event,
        numRegardedPhases (int): number of regarded pick phases

    Returns:
        _type_: _description_
    """
    RefEventPickRemarks = ["_".join(
        [pick.waveform_id.station_code, pick.phase_hint, readExtra(pick)]) for pick in RefEventPicks]
    ComEventPickRemarks = ["_".join(
        [pick.waveform_id.station_code, pick.phase_hint, readExtra(pick)]) for pick in ComEventPicks]
    newPicks = list(set(ComEventPickRemarks) - set(RefEventPickRemarks))
    numRegardedPhases += len(ComEventPickRemarks) - len(newPicks)
    newPicksIndex = [ComEventPickRemarks.index(
        newPick) for newPick in newPicks]
    return newPicksIndex, numRegardedPhases


def ReviewEventPicks(event):
    """Review event picks

    Args:
        event (obspy.event): an obspy event picks, for updated "reference" event,

    Returns:
        obspy.event: _description_
    """
    IgnoredResource_id = []
    for trialPick in event.picks:
        s, p, w = trialPick.waveform_id.station_code, trialPick.phase_hint, readExtra(
            trialPick)
        for pick in event.picks:
            ss, pp, ww = pick.waveform_id.station_code, pick.phase_hint, readExtra(
                pick)
            if (s, p) == (ss, pp) and float(w) > float(ww):
                IgnoredResource_id.append(trialPick.resource_id)
    for i, pick in enumerate(event.picks):
        if pick.resource_id in IgnoredResource_id:
            event.picks.pop(i)
    for i, amplitude in enumerate(event.amplitudes):
        if amplitude.resource_id in IgnoredResource_id:
            event.amplitudes.pop(i)
    for i, arrival in enumerate(event.preferred_origin().arrivals):
        if arrival.pick_id in IgnoredResource_id:
            event.preferred_origin().arrivals.pop(i)
    return event


def updateEvent(RefEvent, ComEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes):
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
    RefEventAmplitudes = RefEvent.amplitudes
    ComEventAmplitudes = ComEvent.amplitudes
    RefEventArrivals = RefEventPreferredOrigin.arrivals
    ComEventArrivals = ComEventPreferredOrigin.arrivals
    newPicksIndex, numRegardedPhases = manageNewPicks(
        RefEventPicks, ComEventPicks, numRegardedPhases)
    for newPickIndex in newPicksIndex:
        newPickResource_id = ComEventPicks[newPickIndex].resource_id
        newArrivals = [
            arrival for arrival in ComEventArrivals if arrival.pick_id == newPickResource_id]
        newAmplitudes = [
            amplitude for amplitude in ComEventAmplitudes if amplitude.pick_id == newPickResource_id]
        RefEventPicks.append(ComEventPicks[newPickIndex])
        RefEventArrivals.extend(newArrivals)
        RefEventAmplitudes.extend(newAmplitudes)
        if "P" in ComEventPicks[newPickIndex].phase_hint:
            numNewPhaseP += 1
        elif "S" in ComEventPicks[newPickIndex].phase_hint:
            numNewPhaseS += 1
        elif "AML" in ComEventPicks[newPickIndex].phase_hint:
            numNewAmplitudes += 1
    RefEvent.picks = RefEventPicks
    RefEvent.preferred_origin().arrivals = RefEventArrivals
    if configs["ReplaceNewPhaseWithHigherWeight"]:
        RefEvent = ReviewEventPicks(RefEvent)
    return RefEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes


# starts application
numCommonEvents = 0
numNewPhaseP = 0
numNewPhaseS = 0
numNewAmplitudes = 0
numRegardedPhases = 0

print("+++ Starts merging catalogs...")
for eventID, RefEvent in enumerate(tqdm(RefCat)):
    for ComEvent in ComCat:
        if computeDiff(RefEvent, ComEvent):
            RefEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes = updateEvent(
                RefEvent, ComEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes)
            ComEventsCat.append(RefEvent)
            numCommonEvents += 1
            break
    TarCat.append(RefEvent)

TarCat.write("updatedCatalog.dat", format="NORDIC")
if configs["OutputCommonEventsCatalog"]:
    ComEventsCat.write("commonEventsCatalog.dat", format="NORDIC")

# summary
print(
    "+++ Number of common events: {numCommonEvents}".format(numCommonEvents=numCommonEvents))
print(
    "+++ Number of new P phases: {numNewPhaseP}".format(numNewPhaseP=numNewPhaseP))
print(
    "+++ Number of new S phases: {numNewPhaseS}".format(numNewPhaseS=numNewPhaseS))
print("+++ Number of new Amplitudes: {numNewAmplitudes}".format(
    numNewAmplitudes=numNewAmplitudes))
print("+++ Number of regarded phases: {numRegardedPhases}".format(
    numRegardedPhases=numRegardedPhases))
