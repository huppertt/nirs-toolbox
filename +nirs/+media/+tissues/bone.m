function prop = bone( lambda )
  
    mua = [650.0000    0.0398
        651.0000    0.0394
        652.0000    0.0390
        653.0000    0.0386
        654.0000    0.0381
        655.0000    0.0377
        656.0000    0.0373
        657.0000    0.0368
        658.0000    0.0363
        659.0000    0.0358
        660.0000    0.0353
        661.0000    0.0348
        662.0000    0.0343
        663.0000    0.0339
        664.0000    0.0334
        665.0000    0.0329
        666.0000    0.0325
        667.0000    0.0321
        668.0000    0.0317
        669.0000    0.0313
        670.0000    0.0309
        671.0000    0.0306
        672.0000    0.0303
        673.0000    0.0300
        674.0000    0.0296
        675.0000    0.0293
        676.0000    0.0290
        677.0000    0.0287
        678.0000    0.0285
        679.0000    0.0282
        680.0000    0.0279
        681.0000    0.0276
        682.0000    0.0274
        683.0000    0.0271
        684.0000    0.0269
        685.0000    0.0266
        686.0000    0.0264
        687.0000    0.0262
        688.0000    0.0260
        689.0000    0.0257
        690.0000    0.0255
        691.0000    0.0254
        692.0000    0.0252
        693.0000    0.0250
        694.0000    0.0249
        695.0000    0.0247
        696.0000    0.0246
        697.0000    0.0245
        698.0000    0.0244
        699.0000    0.0243
        700.0000    0.0242
        701.0000    0.0241
        702.0000    0.0240
        703.0000    0.0239
        704.0000    0.0238
        705.0000    0.0237
        706.0000    0.0236
        707.0000    0.0235
        708.0000    0.0235
        709.0000    0.0234
        710.0000    0.0233
        711.0000    0.0233
        712.0000    0.0232
        713.0000    0.0232
        714.0000    0.0232
        715.0000    0.0232
        716.0000    0.0231
        717.0000    0.0231
        718.0000    0.0231
        719.0000    0.0232
        720.0000    0.0232
        721.0000    0.0232
        722.0000    0.0232
        723.0000    0.0232
        724.0000    0.0232
        725.0000    0.0233
        726.0000    0.0233
        727.0000    0.0233
        728.0000    0.0234
        729.0000    0.0234
        730.0000    0.0235
        731.0000    0.0235
        732.0000    0.0236
        733.0000    0.0237
        734.0000    0.0238
        735.0000    0.0239
        736.0000    0.0240
        737.0000    0.0241
        738.0000    0.0242
        739.0000    0.0242
        740.0000    0.0243
        741.0000    0.0244
        742.0000    0.0245
        743.0000    0.0245
        744.0000    0.0246
        745.0000    0.0246
        746.0000    0.0247
        747.0000    0.0247
        748.0000    0.0247
        749.0000    0.0248
        750.0000    0.0248
        751.0000    0.0248
        752.0000    0.0249
        753.0000    0.0249
        754.0000    0.0249
        755.0000    0.0249
        756.0000    0.0250
        757.0000    0.0250
        758.0000    0.0250
        759.0000    0.0250
        760.0000    0.0250
        761.0000    0.0251
        762.0000    0.0251
        763.0000    0.0251
        764.0000    0.0251
        765.0000    0.0251
        766.0000    0.0251
        767.0000    0.0251
        768.0000    0.0251
        769.0000    0.0250
        770.0000    0.0250
        771.0000    0.0249
        772.0000    0.0249
        773.0000    0.0249
        774.0000    0.0248
        775.0000    0.0248
        776.0000    0.0248
        777.0000    0.0247
        778.0000    0.0247
        779.0000    0.0247
        780.0000    0.0247
        781.0000    0.0247
        782.0000    0.0246
        783.0000    0.0246
        784.0000    0.0246
        785.0000    0.0246
        786.0000    0.0246
        787.0000    0.0246
        788.0000    0.0246
        789.0000    0.0246
        790.0000    0.0246
        791.0000    0.0246
        792.0000    0.0246
        793.0000    0.0246
        794.0000    0.0246
        795.0000    0.0246
        796.0000    0.0246
        797.0000    0.0246
        798.0000    0.0246
        799.0000    0.0247
        800.0000    0.0247
        801.0000    0.0247
        802.0000    0.0247
        803.0000    0.0247
        804.0000    0.0247
        805.0000    0.0247
        806.0000    0.0247
        807.0000    0.0247
        808.0000    0.0247
        809.0000    0.0247
        810.0000    0.0247
        811.0000    0.0248
        812.0000    0.0248
        813.0000    0.0248
        814.0000    0.0248
        815.0000    0.0248
        816.0000    0.0248
        817.0000    0.0248
        818.0000    0.0248
        819.0000    0.0248
        820.0000    0.0248
        821.0000    0.0248
        822.0000    0.0249
        823.0000    0.0249
        824.0000    0.0249
        825.0000    0.0250
        826.0000    0.0250
        827.0000    0.0251
        828.0000    0.0251
        829.0000    0.0252
        830.0000    0.0252
        831.0000    0.0253
        832.0000    0.0253
        833.0000    0.0254
        834.0000    0.0254
        835.0000    0.0255
        836.0000    0.0255
        837.0000    0.0256
        838.0000    0.0257
        839.0000    0.0257
        840.0000    0.0258
        841.0000    0.0259
        842.0000    0.0260
        843.0000    0.0261
        844.0000    0.0262
        845.0000    0.0263
        846.0000    0.0264
        847.0000    0.0265
        848.0000    0.0265
        849.0000    0.0266
        850.0000    0.0267
        851.0000    0.0268
        852.0000    0.0269
        853.0000    0.0270
        854.0000    0.0271
        855.0000    0.0272
        856.0000    0.0273
        857.0000    0.0274
        858.0000    0.0275
        859.0000    0.0276
        860.0000    0.0277
        861.0000    0.0278
        862.0000    0.0279
        863.0000    0.0281
        864.0000    0.0282
        865.0000    0.0283
        866.0000    0.0284
        867.0000    0.0286
        868.0000    0.0287
        869.0000    0.0288
        870.0000    0.0290
        871.0000    0.0291
        872.0000    0.0293
        873.0000    0.0294
        874.0000    0.0296
        875.0000    0.0298
        876.0000    0.0300
        877.0000    0.0302
        878.0000    0.0303
        879.0000    0.0305
        880.0000    0.0307
        881.0000    0.0310
        882.0000    0.0312
        883.0000    0.0314
        884.0000    0.0316
        885.0000    0.0318
        886.0000    0.0320
        887.0000    0.0323
        888.0000    0.0325
        889.0000    0.0327
        890.0000    0.0330
        891.0000    0.0332
        892.0000    0.0335
        893.0000    0.0338
        894.0000    0.0341
        895.0000    0.0344
        896.0000    0.0347
        897.0000    0.0350
        898.0000    0.0353
        899.0000    0.0356
        900.0000    0.0360
        901.0000    0.0363
        902.0000    0.0366
        903.0000    0.0369
        904.0000    0.0372
        905.0000    0.0374
        906.0000    0.0377
        907.0000    0.0380
        908.0000    0.0382
        909.0000    0.0384
        910.0000    0.0386
        911.0000    0.0389
        912.0000    0.0391
        913.0000    0.0393
        914.0000    0.0395
        915.0000    0.0397
        916.0000    0.0399
        917.0000    0.0401
        918.0000    0.0403
        919.0000    0.0405
        920.0000    0.0407
        921.0000    0.0409
        922.0000    0.0411
        923.0000    0.0413
        924.0000    0.0416
        925.0000    0.0418
        926.0000    0.0421
        927.0000    0.0423
        928.0000    0.0426
        929.0000    0.0429
        930.0000    0.0431
        931.0000    0.0434
        932.0000    0.0437
        933.0000    0.0439
        934.0000    0.0441
        935.0000    0.0444
        936.0000    0.0446
        937.0000    0.0449
        938.0000    0.0451
        939.0000    0.0454
        940.0000    0.0457
        941.0000    0.0460
        942.0000    0.0463
        943.0000    0.0466
        944.0000    0.0470
        945.0000    0.0473
        946.0000    0.0477
        947.0000    0.0480
        948.0000    0.0484
        949.0000    0.0487
        950.0000    0.0491];

    mus = [650.0000    2.7174
        651.0000    2.7078
        652.0000    2.6982
        653.0000    2.6886
        654.0000    2.6790
        655.0000    2.6694
        656.0000    2.6598
        657.0000    2.6502
        658.0000    2.6406
        659.0000    2.6311
        660.0000    2.6216
        661.0000    2.6121
        662.0000    2.6027
        663.0000    2.5932
        664.0000    2.5838
        665.0000    2.5744
        666.0000    2.5650
        667.0000    2.5555
        668.0000    2.5461
        669.0000    2.5367
        670.0000    2.5273
        671.0000    2.5179
        672.0000    2.5085
        673.0000    2.4991
        674.0000    2.4898
        675.0000    2.4805
        676.0000    2.4712
        677.0000    2.4619
        678.0000    2.4527
        679.0000    2.4436
        680.0000    2.4345
        681.0000    2.4255
        682.0000    2.4167
        683.0000    2.4079
        684.0000    2.3993
        685.0000    2.3908
        686.0000    2.3824
        687.0000    2.3743
        688.0000    2.3662
        689.0000    2.3584
        690.0000    2.3508
        691.0000    2.3433
        692.0000    2.3360
        693.0000    2.3289
        694.0000    2.3220
        695.0000    2.3153
        696.0000    2.3088
        697.0000    2.3024
        698.0000    2.2962
        699.0000    2.2902
        700.0000    2.2843
        701.0000    2.2786
        702.0000    2.2730
        703.0000    2.2675
        704.0000    2.2621
        705.0000    2.2567
        706.0000    2.2514
        707.0000    2.2461
        708.0000    2.2409
        709.0000    2.2357
        710.0000    2.2304
        711.0000    2.2252
        712.0000    2.2199
        713.0000    2.2146
        714.0000    2.2092
        715.0000    2.2038
        716.0000    2.1983
        717.0000    2.1928
        718.0000    2.1873
        719.0000    2.1817
        720.0000    2.1761
        721.0000    2.1704
        722.0000    2.1648
        723.0000    2.1591
        724.0000    2.1535
        725.0000    2.1478
        726.0000    2.1422
        727.0000    2.1367
        728.0000    2.1311
        729.0000    2.1256
        730.0000    2.1201
        731.0000    2.1146
        732.0000    2.1092
        733.0000    2.1037
        734.0000    2.0984
        735.0000    2.0930
        736.0000    2.0878
        737.0000    2.0825
        738.0000    2.0774
        739.0000    2.0723
        740.0000    2.0673
        741.0000    2.0623
        742.0000    2.0574
        743.0000    2.0526
        744.0000    2.0479
        745.0000    2.0433
        746.0000    2.0387
        747.0000    2.0342
        748.0000    2.0297
        749.0000    2.0254
        750.0000    2.0211
        751.0000    2.0169
        752.0000    2.0127
        753.0000    2.0087
        754.0000    2.0046
        755.0000    2.0007
        756.0000    1.9967
        757.0000    1.9929
        758.0000    1.9890
        759.0000    1.9852
        760.0000    1.9814
        761.0000    1.9776
        762.0000    1.9738
        763.0000    1.9700
        764.0000    1.9662
        765.0000    1.9624
        766.0000    1.9585
        767.0000    1.9547
        768.0000    1.9508
        769.0000    1.9469
        770.0000    1.9431
        771.0000    1.9391
        772.0000    1.9352
        773.0000    1.9313
        774.0000    1.9274
        775.0000    1.9235
        776.0000    1.9197
        777.0000    1.9158
        778.0000    1.9120
        779.0000    1.9083
        780.0000    1.9046
        781.0000    1.9010
        782.0000    1.8974
        783.0000    1.8939
        784.0000    1.8905
        785.0000    1.8871
        786.0000    1.8837
        787.0000    1.8805
        788.0000    1.8772
        789.0000    1.8741
        790.0000    1.8709
        791.0000    1.8678
        792.0000    1.8648
        793.0000    1.8617
        794.0000    1.8587
        795.0000    1.8557
        796.0000    1.8528
        797.0000    1.8498
        798.0000    1.8468
        799.0000    1.8439
        800.0000    1.8409
        801.0000    1.8379
        802.0000    1.8350
        803.0000    1.8320
        804.0000    1.8290
        805.0000    1.8261
        806.0000    1.8231
        807.0000    1.8201
        808.0000    1.8171
        809.0000    1.8141
        810.0000    1.8111
        811.0000    1.8080
        812.0000    1.8050
        813.0000    1.8020
        814.0000    1.7989
        815.0000    1.7959
        816.0000    1.7929
        817.0000    1.7898
        818.0000    1.7868
        819.0000    1.7838
        820.0000    1.7808
        821.0000    1.7778
        822.0000    1.7748
        823.0000    1.7719
        824.0000    1.7690
        825.0000    1.7660
        826.0000    1.7631
        827.0000    1.7602
        828.0000    1.7573
        829.0000    1.7543
        830.0000    1.7514
        831.0000    1.7484
        832.0000    1.7454
        833.0000    1.7424
        834.0000    1.7394
        835.0000    1.7363
        836.0000    1.7332
        837.0000    1.7301
        838.0000    1.7270
        839.0000    1.7238
        840.0000    1.7206
        841.0000    1.7174
        842.0000    1.7142
        843.0000    1.7110
        844.0000    1.7077
        845.0000    1.7045
        846.0000    1.7013
        847.0000    1.6981
        848.0000    1.6949
        849.0000    1.6917
        850.0000    1.6885
        851.0000    1.6854
        852.0000    1.6823
        853.0000    1.6793
        854.0000    1.6763
        855.0000    1.6733
        856.0000    1.6704
        857.0000    1.6676
        858.0000    1.6647
        859.0000    1.6619
        860.0000    1.6592
        861.0000    1.6565
        862.0000    1.6538
        863.0000    1.6512
        864.0000    1.6486
        865.0000    1.6461
        866.0000    1.6436
        867.0000    1.6411
        868.0000    1.6386
        869.0000    1.6362
        870.0000    1.6337
        871.0000    1.6313
        872.0000    1.6289
        873.0000    1.6266
        874.0000    1.6242
        875.0000    1.6219
        876.0000    1.6195
        877.0000    1.6171
        878.0000    1.6147
        879.0000    1.6123
        880.0000    1.6098
        881.0000    1.6072
        882.0000    1.6045
        883.0000    1.6018
        884.0000    1.5989
        885.0000    1.5959
        886.0000    1.5929
        887.0000    1.5898
        888.0000    1.5865
        889.0000    1.5832
        890.0000    1.5798
        891.0000    1.5763
        892.0000    1.5728
        893.0000    1.5691
        894.0000    1.5654
        895.0000    1.5617
        896.0000    1.5578
        897.0000    1.5539
        898.0000    1.5500
        899.0000    1.5460
        900.0000    1.5419
        901.0000    1.5378
        902.0000    1.5336
        903.0000    1.5294
        904.0000    1.5252
        905.0000    1.5209
        906.0000    1.5167
        907.0000    1.5124
        908.0000    1.5082
        909.0000    1.5039
        910.0000    1.4997
        911.0000    1.4954
        912.0000    1.4911
        913.0000    1.4868
        914.0000    1.4824
        915.0000    1.4780
        916.0000    1.4737
        917.0000    1.4692
        918.0000    1.4648
        919.0000    1.4603
        920.0000    1.4559
        921.0000    1.4514
        922.0000    1.4469
        923.0000    1.4424
        924.0000    1.4380
        925.0000    1.4335
        926.0000    1.4291
        927.0000    1.4247
        928.0000    1.4203
        929.0000    1.4159
        930.0000    1.4116
        931.0000    1.4072
        932.0000    1.4029
        933.0000    1.3986
        934.0000    1.3943
        935.0000    1.3900
        936.0000    1.3857
        937.0000    1.3815
        938.0000    1.3773
        939.0000    1.3731
        940.0000    1.3689
        941.0000    1.3647
        942.0000    1.3605
        943.0000    1.3563
        944.0000    1.3521
        945.0000    1.3479
        946.0000    1.3437
        947.0000    1.3395
        948.0000    1.3353
        949.0000    1.3312
        950.0000    1.3270];

    mua = interp1(mua(:,1),mua(:,2),lambda);
    mus = interp1(mus(:,1),mus(:,2),lambda);
    ri = 1.45;

    prop = nirs.media.OpticalProp(mua,mus,lambda,ri);

end

