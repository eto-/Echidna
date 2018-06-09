void SetAliases() {
     bxtree->SetAlias ("rhits","laben.raw_hits");
     bxtree->SetAlias ("nrhits","laben.n_raw_hits");
     bxtree->SetAlias ("dhits","laben.decoded_hits");
     bxtree->SetAlias ("ndhits","laben.n_decoded_hits");
     bxtree->SetAlias ("dhits.ttime","laben.decoded_hits.raw_time - laben.trigger_time");
     bxtree->SetAlias ("dhits.ltime","laben.decoded_hits.raw_time - laben.laser_time");
     bxtree->SetAlias ("dhits.ctime1","laben.decoded_hits.raw_time - laben.clusters[0].start_time");
     bxtree->SetAlias ("dhits.ctime2","laben.decoded_hits.raw_time - laben.clusters[1].start_time");
     bxtree->SetAlias ("dhits.ctime3","laben.decoded_hits.raw_time - laben.clusters[2].start_time");
     bxtree->SetAlias ("chits","laben.clustered_hits");
     bxtree->SetAlias ("nchits","laben.n_clustered_hits");
     bxtree->SetAlias ("cls","laben.clusters");
     bxtree->SetAlias ("ttype","trigger.trgtype");


     bxtree->SetAlias ("r_msk_1","sqrt(laben.clusters[0].position_msk.x**2 + laben.clusters[0].position_msk.y**2 + laben.clusters[0].position_msk.z**2)");
     bxtree->SetAlias ("r_msk_2","sqrt(laben.clusters[1].position_msk.x**2 + laben.clusters[1].position_msk.y**2 + laben.clusters[1].position_msk.z**2)");
     bxtree->SetAlias ("r_msk_3","sqrt(laben.clusters[2].position_msk.x**2 + laben.clusters[2].position_msk.y**2 + laben.clusters[2].position_msk.z**2)");

     bxtree->SetAlias ("r_mi_1","sqrt(laben.clusters[0].position_mi.x**2 + laben.clusters[0].position_mi.y**2 + laben.clusters[0].position_mi.z**2)");
     bxtree->SetAlias ("r_mi_2","sqrt(laben.clusters[1].position_mi.x**2 + laben.clusters[1].position_mi.y**2 + laben.clusters[1].position_mi.z**2)");
     bxtree->SetAlias ("r_mi_3","sqrt(laben.clusters[2].position_mi.x**2 + laben.clusters[2].position_mi.y**2 + laben.clusters[2].position_mi.z**2)");

     bxtree->SetAlias ("r_lngs_1","sqrt(laben.clusters[0].position_lngs.x**2 + laben.clusters[0].position_lngs.y**2 + laben.clusters[0].position_lngs.z**2)");
     bxtree->SetAlias ("r_lngs_2","sqrt(laben.clusters[1].position_lngs.x**2 + laben.clusters[1].position_lngs.y**2 + laben.clusters[1].position_lngs.z**2)");
     bxtree->SetAlias ("r_lngs_3","sqrt(laben.clusters[2].position_lngs.x**2 + laben.clusters[2].position_lngs.y**2 + laben.clusters[2].position_lngs.z**2)");

     bxtree->SetAlias ("r_dbn_1","sqrt(laben.clusters[0].position_dbn.x**2 + laben.clusters[0].position_dbn.y**2 + laben.clusters[0].position_dbn.z**2)");
     bxtree->SetAlias ("r_dbn_2","sqrt(laben.clusters[1].position_dbn.x**2 + laben.clusters[1].position_dbn.y**2 + laben.clusters[1].position_dbn.z**2)");
     bxtree->SetAlias ("r_dbn_3","sqrt(laben.clusters[2].position_dbn.x**2 + laben.clusters[2].position_dbn.y**2 + laben.clusters[2].position_dbn.z**2)");

}
