expect_scalar <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object))

  # 2. Call expect()
  expect(
    is.numeric(act$val),
    sprintf("%s has type %s, not type 'numeric'.", act$lab, typeof(act$val))
  )

  act$n <- length(act$val)
  expect(
    length(act$n)==1,
    sprintf("%s has length %i, not length 1.", act$lab, act$n)
  )

  expect(
    is.finite(act$val),
    sprintf("%s has non-finite value %i.", act$lab, act$val)
  )

  # 3. Invisibly return the value
  invisible(act$val)
}
