plot_AR_regions <- function(Z, masked_set, sl, sr) {
    df <- data.frame(cbind(pnorm(Z), (1:length(Z) %in% masked_set),
                           sl, sr))
    colnames(df) <- c("U", "masked", "sl", "sr")
    df <- df %>% arrange(desc(masked))
    df$rowid <- 1:length(Z)
    plot <- ggplot(df, aes(rowid, U, color=masked)) +
        geom_point() +
        geom_hline(yintercept=0.5, color="black", linetype="dashed") +
        geom_ribbon(aes(ymin=0, ymax=sl), fill="red", alpha=0.2, color=0) +
        geom_ribbon(aes(ymin=sr, ymax=1), fill="red", alpha=0.2, color=0) +
        geom_ribbon(aes(ymin=0.5-sl, ymax=1.5-sr), fill="blue", alpha=0.2, color=0)
    return(plot)
}
