# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 11:12:49 2018

@author: mirwi
"""


from modeller import *
from modeller.automodel import *    

log.verbose()
env = environ()


env.io.atom_files_directory = ['']

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        



        # secondary structure - chain A
        rsr.add(secondary_structure.alpha(self.residue_range('6:A', '13:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('18:A', '24:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('31:A', '39:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('42:A', '69:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('288:A', '296:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('372:A', '377:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('485:A', '491:A')))

        # secondary structure - chain B
        rsr.add(secondary_structure.alpha(self.residue_range('865:B', '871:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('878:B', '882:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('888:B', '896:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('901:B', '927:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('1144:B', '1152:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('1229:B', '1234:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('1342:B', '1349:B')))

        # cross links

        # Chain A
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:31:A']), mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:76:A'], at['CA:83:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:76:A'], at['CA:84:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:534:A'], at['CA:581:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:765:A'], at['CA:827:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:819:A'], at['CA:827:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:677:A'], at['CA:683:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:827:A'], at['CA:845:A']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:786:A'], at['CA:807:A']),  mean=25.0, stdev=0.5))
        # long range
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:284:A'], at['CA:394:A']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:393:A'], at['CA:620:A']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:534:A'], at['CA:579:A']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:613:A'], at['CA:677:A']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:247:A'], at['CA:326:A']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:613:A'], at['CA:765:A']),  mean=50.0, stdev=0.5))


        # Chain B
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1643:B'], at['CA:1664:B']), mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1299:B'], at['CA:1251:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1349:B'], at['CA:1391:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1391:B'], at['CA:1438:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1391:B'], at['CA:1436:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1534:B'], at['CA:1540:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1546:B'], at['CA:1643:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1622:B'], at['CA:1676:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1691:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1685:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1698:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1678:B'], at['CA:1691:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1692:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1691:B'], at['CA:1698:B']),  mean=10.0, stdev=0.5))
        # long range
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1250:B'], at['CA:1477:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1470:B'], at['CA:1534:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1346:B'], at['CA:1391:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1349:B'], at['CA:1439:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1349:B'], at['CA:1436:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1391:B'], at['CA:1552:B']),  mean=36.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:1546:B'], at['CA:1664:B']),  mean=50.0, stdev=0.5))


        # inter-chain
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:885:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:903:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:31:A'], at['CA:937:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:31:A'], at['CA:885:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:247:A'], at['CA:1250:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:447:A'], at['CA:1299:B']),  mean=25.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:447:A'], at['CA:1477:B']),  mean=25.0, stdev=0.5))
		# long range
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:21:A'], at['CA:937:B']),  mean=36.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:76:A'], at['CA:1043:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:447:A'], at['CA:1251:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:459:A'], at['CA:1477:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:613:A'], at['CA:1174:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:393:A'], at['CA:1104:B']),  mean=50.0, stdev=0.5))
        rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at['CA:677:A'], at['CA:1174:B']),  mean=50.0, stdev=0.5))
		#New
	rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1546:B'], at['CA:1664:B']),  mean=26.4, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1696:B'], at['CA:1691:B']),  mean=12.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1672:B'], at['CA:1534:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1674:B'], at['CA:1534:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1683:B'], at['CA:1691:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1682:B'], at['CA:1691:B']),  mean=12.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1688:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1690:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1303:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1686:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1300:B']),  mean=20.0, stdev=0.5))		
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1302:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1691:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1691:B']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:819:A'], at['CA:821:A']),  mean=30.0, stdev=0.5))	
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:819:A'], at['CA:827:A']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:827:A'], at['CA:829:A']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:827:A'], at['CA:845:A']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1251:B'], at['CA:1299:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1349:B'], at['CA:1391:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1391:B'], at['CA:1436:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1391:B'], at['CA:1438:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1534:B'], at['CA:1540:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1546:B'], at['CA:1643:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1664:B'], at['CA:1667:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1685:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1678:B'], at['CA:1685:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1678:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1686:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1691:B'], at['CA:1692:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1688:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1676:B'], at['CA:1685:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1251:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1691:B']),  mean=30.0, stdev=0.5))	
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:441:A']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:444:A']),  mean=20.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:443:A']),  mean=20.0, stdev=0.5))

        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1334:B'], at['CA:1330:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1672:B'], at['CA:1534:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1674:B'], at['CA:1534:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1683:B'], at['CA:1691:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1682:B'], at['CA:1691:B']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1688:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1690:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1303:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1686:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1303:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1300:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1690:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1686:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1302:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1303:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1690:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1691:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1687:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1303:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1686:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1685:B'], at['CA:1690:B']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:1684:B'], at['CA:1691:B']),  mean=19.0, stdev=0.5))
		
	rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:475:A'], at['CA:474:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:813:A'], at['CA:675:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:815:A'], at['CA:675:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:824:A'], at['CA:832:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:823:A'], at['CA:832:A']),  mean=13.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:829:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:831:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:444:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:825:A'], at['CA:827:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:444:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:441:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:831:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:825:A'], at['CA:827:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:443:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:444:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:831:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:832:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:828:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:444:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:825:A'], at['CA:827:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:826:A'], at['CA:831:A']),  mean=19.0, stdev=0.5))
        rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at['CA:825:A'], at['CA:832:A']),  mean=19.0, stdev=0.5))
		
a = MyModel(env,
            alnfile  = 'Pab.ali',     
            knowns   = ('6mzb.pdb'),              
            sequence = 'Pab',
			assess_methods=(assess.DOPE))              
a.starting_model= 1                 
a.ending_model  = 5               
                                    
a.make()                            