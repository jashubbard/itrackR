context('events, behavioral, drift correction')

z <- itrackr(path=itrackR.data('path'))



test_that("extracting events and merging", {

  z <- set_index(z,varnames=c('Block','Trial'), patterns=c('^BLOCK [0-9]*','^TRIAL [0-9]*'), numeric.only=T)

  z <- find_messages(z,varnames=c('STIMONSET','RESPONSE'), patterns=c('STIMONSET','RESPONSE'), timestamp=T)

  z <- add_behdata(z,itrackR.data('beh'),append=F)

  expect_is(z$beh,'data.frame')
  expect_true(all(c('Block','Trial') %in% names(z$header)))
  expect_true(all(c('STIMONSET','RESPONSE') %in% names(z$header)))
  expect_true(all(c('Block','Trial') %in% names(z$beh)))

  expect_is(z$header$Block, 'numeric')
  expect_is(z$header$Trial, 'numeric')
  expect_is(z$header$STIMONSET, 'numeric')
  expect_is(z$header$RESPONSE, 'numeric')

  orig_x <- z$fixations$gavx
  orig_y <- z$fixations$gavy

  z <- drift_correct(z, vars=c('Block'))
  expect_is(z$transform$fixations,'data.frame')
  expect_false(all(z$fixations$gavx==orig_x))
  expect_false(all(z$fixations$gavy==orig_y))

  z <- undrift(z)
  expect_null(z$transform$fixations)
  expect_true(all(z$fixations$gavx==orig_x))
  expect_true(all(z$fixations$gavy==orig_y))


})


