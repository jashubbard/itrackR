context('basics')

z <- itrackr(path=itrackR.data('path'))

test_that("basic tests", {

  expect_is(z,'itrackR')
  expect_equal(length(z$edfs),2)
  expect_equal(length(z$subs),2)
  expect_equal(z$samples,list())

})




