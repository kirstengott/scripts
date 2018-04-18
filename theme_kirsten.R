
theme_kirsten <- function(rot = TRUE, presentation = FALSE, angle = 45, text_size = 12){
  my_theme <- theme(strip.text.x = element_text(face="bold"),
                    strip.text.y = element_text(face="bold"),
                    strip.background = element_rect(fill = 'grey96'))
  if(rot == TRUE){
    my_theme <- my_theme + theme(axis.text.x = element_text(angle = angle, hjust = 1))
  }
  if(presentation == TRUE){
    my_theme <- my_theme + theme(text = element_text(size = 24))
  }
  theme_bw(base_size = text_size) + my_theme
}
