#' @noRd
redundant_columns <- function(){
    redundant_cols <- c("sample_index","Sample.ID","Flags","Sample.Type","Panel",
                        "Filter.Set","Custom.Annotation","Input.DNA.Mass..ng.",
                        "Input.DNA.Q.Ratio","Sample.Primer.Plate","Sample.Primer",
                        "Subject.ID","Date.of.Sample.Collection","Analysis.Name",
                        "Analysis.ID","Number.of.Read.Pairs","Variants.Before.Filtering",
                        "Variants.After.Filtering","Mapped.Reads","On.Target.Rate",
                        "Sequencing.Depth.Median","Sequencing.Depth.5th.Percentile",
                        "Sequencing.Depth.95th.Percentile","Unique.Depth.Median",
                        "Unique.Depth.5th.Percentile","Unique.Depth.95th.Percentile",
                        "Error.Rate","Fragment.Length.Median","Fragment.Length.5th.Percentile",
                        "Fragment.Length.95th.Percentile","Bases.Within.10.fold.Range.of.Median",
                        "Number.of.Read.Pairs.on.Lane","Bases.at.or.Above.Q30",
                        "Aligned.to.PhiX","Percentage.of.panel.greater.than.300x.coverage",
                        "Analysis.Software.Version","Sequencer.Name","Flowcell.ID",
                        "Analysis.Owner","Sequencer.File.Data.Path","Number.of.Samples",
                        "Sequencer.Start.Date","Analysis.Submission.Date",
                        "Analysis.Submission.Time","Analysis.Started.Date",
                        "Analysis.Started.Time","Analysis.Completed.Date",
                        "Analysis.Completed.Time","CSV.Report.Generation.Date",
                        "CSV.Report.Generation.Time","Isolated.DNA.Mass..ng.",
                        "Plasma.Volume..mL.","Sample.Adapter",
                        "Theoretical.Sensitivity.Median",
                        "Theoretical.Sensitivity.5th.Percentile",
                        "Theoretical.Sensitivity.95th.Percentile")
    return(redundant_cols)
}
