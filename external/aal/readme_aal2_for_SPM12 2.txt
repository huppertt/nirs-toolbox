Anatomical Automatic Labeling 2 (AAL2) for SPM12
Date 25/08/2015
E-mail: aal2.gin.brainatlas@gmail.com

This archive provide an alternative parcellation of the orbitofrontal cortex for the automated anatomical labeling atlas of Tzourio-Mazoyer et al. (Tzourio-Mazoyer et al., 2002). The new parcellation of the orbitofrontal cortex follows the description provided by Chiavaras, Petrides, and colleagues (Chiavaras and Petrides, 2000; Chiavaras et al., 2001).

The AAL2 atlas is described in the following publication:
Rolls ET, Joliot M, Tzourio-Mazoyer N (2015) Implementation of a new parcellation of the orbitofrontal cortex in the automated anatomical labeling atlas. NeuroImage 10.1016/j.neuroimage.2015.07.075.


Anatomical Automatic Labeling (AAL2) is a package for the anatomical labeling of functional brain mapping experiments. It is an in-house package made by the Neurofunctional Imaging Group (GIN, UMR5296, Bordeaux, France), which is available to the scientific community as copyright freeware under the terms of the GNU General Public License.


The archive file aal2.tar.gz includes 2 sets of files: 
- The first set is to be used with the AAL software for SPM12 (provided through http://www.gin.cnrs.fr/AAL2)
ROI_MNI_V5_vol.mat
ROI_MNI_V5_List.mat
ROI_MNI_V5_Border.mat
ROI_MNI_V5_vol.txt
ROI_MNI_V5.xml
ROI_MNI_V5.nii
AAL2.xml
AAL2.nii

- The second set is to be used with the mricron (http://www.mccauslandcenter.sc.edu/mricro/mricron) software: 
aal2.nii.gz
aal2.nii.lut
aal2.nii.txt

Note 1: As shown in Table 2 of (Rolls, Joliot and Tzourio-Mazoyer 2015), eight new labels have been provided,
OFCmed_L
OFCmed_R
OFCant_L
OFCant_R
OFCpost_L
OFCpost_R
OFClat_l
OFClat_R

Note 2: As shown in Table 2 of (Rolls, Joliot and Tzourio-Mazoyer 2015) six labels have been modified. This was chosen to emphasis that those regions, already present in the first version of aal (aal and ROI_MNI_V4), have been modified in the aal2 parcellation. 
Frontal_Sup_L changes to Frontal_Sup_2_L
Frontal_Sup_R changes to Frontal_Sup_2_R
Frontal_Mid_L changes to Frontal_Mid_2_L
Frontal_Mid_R changes to Frontal_Mid_2_R
Frontal_Inf_Orb_L changes to Frontal_Inf_Orb_2_L
Frontal_Inf_Orb_R changes to Frontal_Inf_Orb_2_R

Note 3: As shown in Table 2 of (Rolls, Joliot and Tzourio-Mazoyer 2015) six labels have been modified without regional modification:
Cingulum_Ant_L changes to Cingulate_Ant_L
Cingulum_Ant_R changes to Cingulate_Ant_R
Cingulum_Mid_L changes to Cingulate_Mid_L
Cingulum_Mid_R changes to Cingulate_Mid_R
Cingulum_Post_L changes to Cingulate_Post_L
Cingulum_Post_R changes to Cingulate_Post_R

Note 4: aal2 and ROI_MNI_V5 differ only in the index see below.
anatomical-name	index-aal2	index-ROI_MNI_V5
Precentral_L	1	2001
Precentral_R	2	2002
Frontal_Sup_2_L	3	2101
Frontal_Sup_2_R	4	2102
Frontal_Mid_2_L	5	2201
Frontal_Mid_2_R	6	2202
Frontal_Inf_Oper_L	7	2301
Frontal_Inf_Oper_R	8	2302
Frontal_Inf_Tri_L	9	2311
Frontal_Inf_Tri_R	10	2312
Frontal_Inf_Orb_2_L	11	2321
Frontal_Inf_Orb_2_R	12	2322
Rolandic_Oper_L	13	2331
Rolandic_Oper_R	14	2332
Supp_Motor_Area_L	15	2401
Supp_Motor_Area_R	16	2402
Olfactory_L	17	2501
Olfactory_R	18	2502
Frontal_Sup_Medial_L	19	2601
Frontal_Sup_Medial_R	20	2602
Frontal_Med_Orb_L	21	2611
Frontal_Med_Orb_R	22	2612
Rectus_L	23	2701
Rectus_R	24	2702
OFCmed_L	25	2801
OFCmed_R	26	2802
OFCant_L	27	2811
OFCant_R	28	2812
OFCpost_L	29	2821
OFCpost_R	30	2822
OFClat_L	31	2831
OFClat_R	32	2832
Insula_L	33	3001
Insula_R	34	3002
Cingulate_Ant_L	35	4001
Cingulate_Ant_R	36	4002
Cingulate_Mid_L	37	4011
Cingulate_Mid_R	38	4012
Cingulate_Post_L	39	4021
Cingulate_Post_R	40	4022
Hippocampus_L	41	4101
Hippocampus_R	42	4102
ParaHippocampal_L	43	4111
ParaHippocampal_R	44	4112
Amygdala_L	45	4201
Amygdala_R	46	4202
Calcarine_L	47	5001
Calcarine_R	48	5002
Cuneus_L	49	5011
Cuneus_R	50	5012
Lingual_L	51	5021
Lingual_R	52	5022
Occipital_Sup_L	53	5101
Occipital_Sup_R	54	5102
Occipital_Mid_L	55	5201
Occipital_Mid_R	56	5202
Occipital_Inf_L	57	5301
Occipital_Inf_R	58	5302
Fusiform_L	59	5401
Fusiform_R	60	5402
Postcentral_L	61	6001
Postcentral_R	62	6002
Parietal_Sup_L	63	6101
Parietal_Sup_R	64	6102
Parietal_Inf_L	65	6201
Parietal_Inf_R	66	6202
SupraMarginal_L	67	6211
SupraMarginal_R	68	6212
Angular_L	69	6221
Angular_R	70	6222
Precuneus_L	71	6301
Precuneus_R	72	6302
Paracentral_Lobule_L	73	6401
Paracentral_Lobule_R	74	6402
Caudate_L	75	7001
Caudate_R	76	7002
Putamen_L	77	7011
Putamen_R	78	7012
Pallidum_L	79	7021
Pallidum_R	80	7022
Thalamus_L	81	7101
Thalamus_R	82	7102
Heschl_L	83	8101
Heschl_R	84	8102
Temporal_Sup_L	85	8111
Temporal_Sup_R	86	8112
Temporal_Pole_Sup_L	87	8121
Temporal_Pole_Sup_R	88	8122
Temporal_Mid_L	89	8201
Temporal_Mid_R	90	8202
Temporal_Pole_Mid_L	91	8211
Temporal_Pole_Mid_R	92	8212
Temporal_Inf_L	93	8301
Temporal_Inf_R	94	8302
Cerebelum_Crus1_L	95	9001
Cerebelum_Crus1_R	96	9002
Cerebelum_Crus2_L	97	9011
Cerebelum_Crus2_R	98	9012
Cerebelum_3_L	99	9021
Cerebelum_3_R	100	9022
Cerebelum_4_5_L	101	9031
Cerebelum_4_5_R	102	9032
Cerebelum_6_L	103	9041
Cerebelum_6_R	104	9042
Cerebelum_7b_L	105	9051
Cerebelum_7b_R	106	9052
Cerebelum_8_L	107	9061
Cerebelum_8_R	108	9062
Cerebelum_9_L	109	9071
Cerebelum_9_R	110	9072
Cerebelum_10_L	111	9081
Cerebelum_10_R	112	9082
Vermis_1_2	113	9100
Vermis_3	114	9110
Vermis_4_5	115	9120
Vermis_6	116	9130
Vermis_7	117	9140
Vermis_8	118	9150
Vermis_9	119	9160
Vermis_10	120	9170



> How to install the software

1) Copy the archive to the chosen location (e.g. SPM12 installed in /usr/local/soft)
unix> cp aal2_for_SPM12.tar.gz /usr/local/soft/spm12/toolbox
unix> cd /usr/local/soft/spm12/toolbox

2) Gunzip and untar the archive will create an aal directory
unix> tar -zxvf aal2_for_SPM12.tar.gz

3) Add this directory to your Matlab path and copy 2 files in your SPM12/atlas directory
unix> setenv MATLABPATH $MATLABPATH:/usr/local/soft/spm12/toolbox/aal
unix> cp /usr/local/soft/spm12/toolbox/aal/atlas/AAL2.nii /usr/local/soft/spm12/atlas
unix> cp /usr/local/soft/spm12/toolbox/aal/atlas/AAL2.xml /usr/local/soft/spm12/atlas

4) To install AAL2 in mricron software (e.g. mricron installed in /usr/local/soft)
unix> cp /usr/local/soft/spm12/toolbox/aal/aal2.nii.gz /usr/local/soft/mricron/templates
unix> cp /usr/local/soft/spm12/toolbox/aal/aal2.nii.lut /usr/local/soft/mricron/templates
unix> cp /usr/local/soft/spm12/toolbox/aal/aal2.nii.txt /usr/local/soft/mricron/templates


> How to use the software

1) First option: launch AAL from SPM12:
unix> matlab
>> spm fmri
Select the desired contrast, mask, probability and extent threshold like in the regular spm_result.
In the SPM12 Menu window: Toolbox / aal
Then choose a labeling procedure as below from the instructions 3).

2) Second option: launch AAL from the Matlab command window:
>> aal
Select the desired contrast, mask, probability and extent threshold like in the regular spm_result.
Then choose a labeling procedure as below from the instructions 3).

3) Choose a labeling procedure. The 3 choices are explained and documented in
the NeuroImage paper (Tzourio-Mazoyer et al., 2002):

    Local maxima labeling
    Extended local maxima labeling
    Cluster labeling

4) For "Extended local maxima labeling" input the local maxima radius of the
sphere in millimeters (default 10 mm).

5) Select the anatomical parcellation database
In /usr/local/soft/spm12/toolbox/aal
The file: ROI_MNI_V5.img

6) Results
6.1) Local maxima labeling
For each local maxima:
-coordinates in mm x,y,z
-anatomical label (see below)
-distance in millimeter to this region. If the local maxima is inside a region this distance is null (0.00). If the local maxima is outside the parcellation the nearest region name is displayed in the previous column and the shortest distance from the local maxima to this region is listed (exp: 2.30 mm)
-anatomical label of the local maxima to the second nearest region
-shortest distance of the local maxima to the second nearest region
-anatomical label of the local maxima to the third nearest region
-shortest distance of the local maxima to the third nearest region

6.2) Extended local maxima labeling
Each local maxima is supposed to be a 10 mm (if the default is used) spherical region. The intersection of this volume and the AVOI is computed and the result sorted in a descending order according the percentage of overlap (exp: a result of Postcentral_L 100 % indicates that the 10mm radius region surrounding the local maxima is fully included in the Postcentral_L region)
For each local maximum:
-coordinates in mm x,y,z
-list of anatomical label and percentage of overlap. Percentage less than 1% is not listed. If part of the region is outside the parcellation the anatomical label will list "OUTSIDE".

6.3) Cluster labeling
The intersection of each cluster and the AVOI is computed and the result sorted in a descending order according the percentage of overlap.
For each local maximum:
-coordinates in mm x,y,z of the most significant local maximum of the cluster
-list of anatomical label and percentage of overlap. Percentage less than 1% is not listed. If part of the region is outside the parcellation the anatomical label will list "OUTSIDE".
 
Example: a result of
-44 -22 - 56       Postcentral_L     55.00
          		Precentral_L     31.00
          		OUTSIDE         10.00
          		Parietal_Sup_L      5.00
indicates that:

55% of the cluster volume is included in the Postcentral_L region
31% of the cluster volume is included in the Precentral_L region
10% of the cluster volume is outside the parcellation
5% of the cluster volume is included in the Parietal_Sup_L region

_

There is also a third option to get the label:
unix> matlab
>> spm fmri
Select the desired contrast, mask, probability and extent threshold like in the regular spm_result.
In the SPM12 Results window: Atlas / Label using / AAL2
Then you get the label with a right click on the coordinates in the Graphic window.



Bibliography
Chiavaras MM, Petrides M (2000) Orbitofrontal sulci of the human and macaque monkey brain. J Comp Neurol 422:35-54.
Chiavaras MM, LeGoualher G, Evans A, Petrides M (2001) Three-dimensional probabilistic atlas of the human orbitofrontal sulci in standardized stereotaxic space. Neuroimage 13:479-496.
Rolls ET, Joliot M, Tzourio-Mazoyer N (2015) Implementation of a new parcellation of the orbitofrontal cortex in the automated anatomical labeling atlas. Neuroimage 10.1016/j.neuroimage.2015.07.075.
Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O, Delcroix N, Mazoyer B, Joliot M (2002) Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI single-subject brain. Neuroimage 15:273-289.



History
Version vbeta0 Date 01/01/2000
Version vbeta1 Date 01/01/2001
Version vbeta2 Date 26/02/2002
Version vbeta1 Date 17/09/2003
Version v1 Date 14/06/2010
Version v4 Date 25/08/2015
Version v5 Date 25/08/2015
