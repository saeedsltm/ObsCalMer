import warnings
from json import load

import obspy as obs
from obspy.core.event import Catalog
from obspy.geodetics.base import gps2dist_azimuth as gps
from tqdm import tqdm

warnings.filterwarnings("ignore")

class ObspyCatalogMerger():
    def __init__(self):
        with open("config.json") as f:
            self.configs = load(f)        
        self.RefCat = obs.read_events(self.configs["ReferenceCatalogFileName"])
        self.ComCat = obs.read_events(self.configs["ComparedCatalogFileName"])
        self.TarCat = Catalog()
        self.ComEventsCat = Catalog()
        self.OriginTimeShift = self.configs["OriginTimeShift"]
        self.EpicShift = self.configs["EpicentralShift"]

    def ComputeDiff(self, RefEvent, ComEvent):
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
        if dT < self.OriginTimeShift and dD < self.EpicShift:
            return True
        else:
            return False


    def ReadExtra(self, pick):
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

    def FilterPhase(self, pick):
        """Filter phase based on user request

        Args:
            pick (event.pick): an obspy event pick

        Returns:
            bool: true if user decides P,S or Amplitude to be added
        """
        if not self.configs["AddPhaseP"] and "P" in pick.phase_hint:
            return False
        if not self.configs["AddPhaseS"] and "S" in pick.phase_hint:
            return False
        if not self.configs["AddAmplitude"] and "AML" in pick.phase_hint:
            return False            
        return True

    def ManageNewPicks(self, RefEventPicks, ComEventPicks, numRegardedPhases):
        """Manage new picks

        Args:
            RefEventPicks (event.picks): an obspy event picks, for "reference" event,
            ComEventPicks (event.picks): an obspy event picks, for "compared" event,
            numRegardedPhases (int): number of regarded pick phases

        Returns:
            _type_: _description_
        """
        RefEventPickRemarks = ["_".join(
            [pick.waveform_id.station_code, pick.phase_hint, self.ReadExtra(pick)]) for pick in RefEventPicks]
        ComEventPickRemarks = ["_".join(
            [pick.waveform_id.station_code, pick.phase_hint, self.ReadExtra(pick)]) if self.FilterPhase(pick) else None for pick in ComEventPicks]
        newPicks = list(filter(None, set(ComEventPickRemarks) - set(RefEventPickRemarks)))
        numRegardedPhases += (len(ComEventPickRemarks) - len(newPicks))
        newPicksIndex = [ComEventPickRemarks.index(
            newPick) for newPick in newPicks]
        return newPicksIndex, numRegardedPhases


    def ReviewEventPicks(self, Event):
        """Review event picks

        Args:
            Event (obspy.event): an obspy event picks, for updated "reference" event,

        Returns:
            obspy.event: _description_
        """
        IgnoredResource_id = []
        for trialPick in Event.picks:
            s, p, w = trialPick.waveform_id.station_code, trialPick.phase_hint, self.ReadExtra(
                trialPick)
            for pick in Event.picks:
                ss, pp, ww = pick.waveform_id.station_code, pick.phase_hint, self.ReadExtra(
                    pick)
                if (s, p) == (ss, pp) and float(w) > float(ww):
                    IgnoredResource_id.append(trialPick.resource_id)
        for i, pick in enumerate(Event.picks):
            if pick.resource_id in IgnoredResource_id:
                Event.picks.pop(i)
        for i, amplitude in enumerate(Event.amplitudes):
            if amplitude.resource_id in IgnoredResource_id:
                Event.amplitudes.pop(i)
        for i, arrival in enumerate(Event.preferred_origin().arrivals):
            if arrival.pick_id in IgnoredResource_id:
                Event.preferred_origin().arrivals.pop(i)
        return Event


    def UpdateEvent(self, RefEvent, ComEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes):
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
        newPicksIndex, numRegardedPhases = self.ManageNewPicks(
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
        if self.configs["ReplaceNewPhaseOnlyWithHigherWeight"]:
            RefEvent = self.ReviewEventPicks(RefEvent)
        return RefEvent, numNewPhaseP, numNewPhaseS, numRegardedPhases, numNewAmplitudes

    def MergeCatalog(self):
        """Merge catalogs
        """
        self.numCommonEvents = 0
        self.numNewPhaseP = 0
        self.numNewPhaseS = 0
        self.numNewAmplitudes = 0
        self.numRegardedPhases = 0
        print("+++ Starts merging catalogs...")
        for _, RefEvent in enumerate(tqdm(self.RefCat)):
            for ComEvent in self.ComCat:
                if self.ComputeDiff(RefEvent, ComEvent):
                    RefEvent, self.numNewPhaseP, self.numNewPhaseS, self.numRegardedPhases, self.numNewAmplitudes = self.UpdateEvent(
                        RefEvent, ComEvent, self.numNewPhaseP, self.numNewPhaseS, self.numRegardedPhases, self.numNewAmplitudes)
                    self.ComEventsCat.append(RefEvent)
                    self.numCommonEvents += 1
                    break
            self.TarCat.append(RefEvent)

        self.TarCat.write("updatedCatalog.dat", format="NORDIC")
        if self.configs["OutputCommonEventsCatalog"]:
            self.ComEventsCat.write("commonEventsCatalog.dat", format="NORDIC")

    def WriteSummary(self):
        """Write a summary of what has been merged
        """
        with open("summary.dat", "w") as f:
            f.write(
                "+++ Number of common events: {numCommonEvents}\n".format(numCommonEvents=self.numCommonEvents))
            f.write(
                "+++ Number of new P phases: {numNewPhaseP}\n".format(numNewPhaseP=self.numNewPhaseP))
            f.write(
                "+++ Number of new S phases: {numNewPhaseS}\n".format(numNewPhaseS=self.numNewPhaseS))
            f.write("+++ Number of new amplitudes: {numNewAmplitudes}\n".format(
                numNewAmplitudes=self.numNewAmplitudes))
            f.write("+++ Number of regarded phases: {numRegardedPhases}\n".format(
                numRegardedPhases=self.numRegardedPhases))
            f.write("+++ updated catalog file: 'updatedCatalog.dat'\n")
            f.write("+++ common-events-only catalog file: 'commonEventsCatalog.dat'\n")
            print("+++ summary has been written to file 'summary.dat'.")

if __name__ == "__main__":
    app = ObspyCatalogMerger()
    app.MergeCatalog()
    app.WriteSummary()