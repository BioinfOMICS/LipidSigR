# Test for error input

    Code
      char_2wayAnova(processed_se = processed_se_multiGroup, ratio_transform = NULL,
        char_transform = "log10")
    Condition
      Error in `convert_sp2ratio()`:
      ! The 'transform' parameter must be either 'none' or 'log2'.

---

    Code
      char_2wayAnova(processed_se = processed_se_multiGroup, ratio_transform = "log2",
        char_transform = "log2")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      char_2wayAnova(processed_se = sub, ratio_transform = "log2", char_transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `char_2wayAnova()`:
      ! This function requires each group to have at least two samples.

