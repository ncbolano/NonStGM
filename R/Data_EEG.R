csv_path = "C:/Users/noahb/Downloads/EEG_full.csv"
EEG_full = read.csv(
  csv_path,
  header = TRUE,
  check.names = FALSE,
  na.strings = c("NA", ""),
  colClasses = "numeric"
)
expected_columns = c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2", "SP1(CORRUPTED)", "SP2(CORRUPTED)", "F7", "F8", "T3", "T4", "T5", "T6", "Fz", "Cz(REFERENCE)", "Pz", "EKG", "EKG2", "LOC", "ROC", "CHN", "CHN2", "ch28", "ch29", "ch30", "ch31")
