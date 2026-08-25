# lipid_profiling can handle invalid char

    Code
      lipid_profiling(processed_se, char = "lipid")
    Condition
      Error in `lipid_profiling()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      lipid_profiling(processed_se, char = NULL)
    Condition
      Error in `lipid_profiling()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

