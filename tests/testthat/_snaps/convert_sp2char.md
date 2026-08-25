# convert_sp2char return error message for special condition or incorrect input

    Code
      convert_sp2char(processed_se, transform = "log20")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      convert_sp2char(sub, transform = "log10")
    Condition
      Error in `convert_sp2char()`:
      ! No enough data.

