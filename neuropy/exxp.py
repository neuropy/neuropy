"""Parse stimulus trial information from Expo .xml file. This is a work in progress..."""

import xml.etree.ElementTree as ET

## might need to add the </ExpoXData> close tag at end of exxp xml file for parsing to succeed

fname = '/home/mspacek/data/blab/ephys/PVCre_0108/s02/PVCre_0108_s02_20160401_06.xml'

tree = ET.parse(fname)
root = tree.getroot()
assert root.tag == 'ExpoXData'

blocks = root.find('Blocks')
slots = root.find('Slots')
environment = root.find('Environment')
passes = root.find('Passes')
npasses = int(passes.get('NumOfPasses'))
assert len(passes) == npasses

passkeys = ['ID', 'SlotID', 'BlockID', 'StartTime', 'EndTime']

trials = np.zeros((npasses, 5), dtype=np.int64) # slotid, blockid, tstart, tend, dt
for i, pas in enumerate(passes):
    passi, sloti, blocki, tstart, tend = np.int64([ pas.attrib[key] for key in passkeys ])
    assert passi == i # make sure no passes are missing in the .xml file
    # convert timestamps from units of 100 us to us:
    tstart *= 100
    tend *= 100
    trials[i] = sloti, blocki, tstart, tend, tend-tstart

# make sure all timestamps are in increasing order:
assert (np.diff(trials[:, 2]) > 0).all()
assert (np.diff(trials[:, 3]) > 0).all()
