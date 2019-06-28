context("Parser functions")
test_that("Nexus file can be parsed", {
  filename <- 'test-parse-nexus.nexus'
  read <- ReadCharacters(filename)
  expect_equal(192, ncol(read))
  expect_equal(80, nrow(read))
  expect_equal("Wiwaxia", rownames(read)[4])
  expect_equal("(01)", as.character(read[1, 27]))
})

test_that("TNT trees parsed correctly", {
  trees <- ReadTntTree('test-tnt-tree.tre')
  expect_equal(2, length(trees))
  expect_equal(32, ConsensusWithout(trees, 'Paterimitra')$Nnode)
  
  fromLabels <- ReadTntTree('test-tnt-tree.tre', tipLabels = c('Dailyatia', 'Novocrania', 'Craniops', 'Ussunia', 'Gasconsia', 'Heliomedusa_orienta', 'Micrina', 'Mickwitzia_muralensis', 'Micromitra', 'Askepasma_toddense', 'Pelagodiscus_atlanticus', 'Lingula', 'Eoobolus', 'Clupeafumosus_socialis', 'Phoronis', 'Eccentrotheca', 'Yuganotheca_elegans', 'Longtancunella_chengjiangensis', 'Paterimitra', 'Lingulellotreta_malongensis', 'Acanthotretella', 'Lingulosacculus', 'Pedunculotheca_diania', 'Haplophrentis_carinatus', 'Tomteluva_perturbata', 'Salanygolina', 'Mummpikia_nuda', 'Alisina', 'Coolinia_pecten', 'Antigonambonites_planus', 'Kutorgina_chengjiangensis', 'Nisusia_sulcata', 'Glyptoria', 'Orthis', 'Terebratulina'))
  expect_identical(trees, fromLabels)
  
  namedLabels <- ReadTntTree('test-tnt-namedtree.tre')[[1]]$tip.label
  expect_equal('Flustra', namedLabels[1])
  expect_equal(74L, length(namedLabels))
})
