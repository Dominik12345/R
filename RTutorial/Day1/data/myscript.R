setClass(Class = "Song" , 
         representation(title = 'character',
                        composer = 'character',
                        genre = 'character',
                        year = 'integer',
                        artist = 'character',
                        link = 'character'),
         prototype(title = NA_character_,
                   composer = NA_character_,
                   genre = NA_character_,
                   year = NA_integer_,
                   artist = NA_character_,
                   link = NA_character_) 
)

play <- function(object) {
  UseMethod('play')
}

setMethod('play',signature = 'Song' ,
          function(object) {
            browseURL(object@link)
          }
)

# Test
song1 <- new('Song', title = 'Let It Be', composer = 'Paul McCartney', 
             genre = 'Rock Ballade' , year = as.integer(1970), 
             artist = 'The Beatles', 
             link = 'https://www.youtube.com/watch?v=2xDzVZcqtYI')
song2 <- new('Song' , title = 'Highway to Hell', 
              composer = 'Angus Young, Malcolm Young and Bon Scott', 
              genre = 'Hard Rock' , year = as.integer(1979), 
              artist = 'AC/DC',
              link = 'https://www.youtube.com/watch?v=l482T0yNkeo')
