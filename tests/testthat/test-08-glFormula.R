context("extract for trees")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

glCall <- quote(
  stan4bart:::glFormula(
    formula = formula, data = df,
     control = stan4bart:::make_glmerControl(),
     na.action = getOption("na.action", "na.omit")
   )
)

test_that("glFormula correctly identifies component variables", {

  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X3 - X4 - g.1 - g.2 + X3:X4 + (1 + X4 | g.1) + (1 | g.2))
  
  expect_warning(glResult <- eval(glCall))
  

  expect_setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2"))
  expect_setequal(colnames(glResult$bartData@x), c("X1", "X2"))

  expect_setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 1:10), "z"))
  expect_setequal(colnames(glResult$X), c("X1", "X2", "X3:X4", paste0("X", 5:10), "z"))

  expect_setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4"))
  expect_setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1"))
})

test_that("glFromula removes offset from bart frame", {
  offset <- rep(c(1:10), length.out = nrow(df))

  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X1 - X2 - g.1 - g.2 + (1 + X4 | g.1) + (1 | g.2))
  glCall$offset <- quote(offset)
  
  glResult <- eval(glCall)

  expect_setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2"))
  expect_setequal(colnames(glResult$bartData@x), c("X1", "X2"))

  expect_setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 3:10), "z"))
  expect_setequal(colnames(glResult$X), c(paste0("X", 3:10), "z"))

  expect_setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4"))
  expect_setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1"))

  df$offset <- offset
  rm(offset)
  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X1 - X2 - g.1 - g.2 - offset + (1 + X4 | g.1) + (1 | g.2))
  glResult <- eval(glCall)

  expect_setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2"))
  expect_setequal(colnames(glResult$bartData@x), c("X1", "X2"))

  expect_setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 3:10), "z"))
  expect_setequal(colnames(glResult$X), c(paste0("X", 3:10), "z"))

  expect_setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4"))
  expect_setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1"))
})

test_that("works when parametric component contains only nested factors", {
  glCall$formula <- quote(y ~ bart(. - g.1 - g.2 - X4 - z) + (1 | g.1:g.2))
  glResult <- eval(glCall)

  expect_setequal(attr(glResult$terms, "varnames.bart"), paste0("X", c(1:3, 5:10)))
  expect_setequal(colnames(glResult$bartData@x), paste0("X", c(1:3, 5:10)))
  
  expect_equal(attr(glResult$terms, "varnames.fixed"), character())
  expect_equal(dim(glResult$X), c(nrow(df), 0L))

  expect_setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2"))
  expect_equal(names(glResult$reTrms$Ztlist), "1 | g.1:g.2")
})
