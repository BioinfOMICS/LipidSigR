# Test that message return for invalid input

    Code
      boxPlot_feature_twoGroup(processed_se_twoGroup, feature = "ABC", ref_group = "ctrl",
        test = "t-test", transform = "log10")
    Condition
      Error in `boxPlot_feature_twoGroup()`:
      ! 'feature' must be one of the lipids in your lipid abundance data; please revise your input.

---

    Code
      boxPlot_feature_twoGroup(processed_se_twoGroup, feature = "Cer 38:1;2",
        ref_group = "ctrl1", test = "t-test", transform = "log10")
    Condition
      Error in `boxPlot_feature_twoGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      boxPlot_feature_twoGroup(processed_se_twoGroup, feature = "Cer 38:1;2",
        ref_group = "ctrl", test = "t test", transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "t-test", "Wilcoxon test"

---

    Code
      boxPlot_feature_twoGroup(processed_se_twoGroup, feature = "Cer 38:1;2",
        ref_group = "ctrl", test = "t-test", transform = "log2")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      boxPlot_feature_twoGroup(processed_se = processed_se_multiGroup, feature = "PE O- 17:0;0_20:3;0",
        ref_group = "ctrl", test = "t-test", transform = "log10")
    Condition
      Error in `boxPlot_feature_twoGroup()`:
      ! This function is used for two groups. If your data consists of multiple groups, please use "boxPlot_feature_multiGroup" for analysis.

