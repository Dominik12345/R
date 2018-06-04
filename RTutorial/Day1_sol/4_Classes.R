###############################################################################
#                                                                             #
# R supports three object structures: S3, S4 and R5                           # 
#                                                                             #  
###############################################################################

# clear the Global Environment
rm(list=ls())

#############################
# 4.1 ad-hoc S3             #
#############################

song1 <- list('title' = 'Let It Be', 'composer' = 'Paul McCartney', 
              'genre' = 'Rock Ballade' , 'year' = 1970, 
              'artist' = 'The Beatles', 
              'link' = 'https://www.youtube.com/watch?v=2xDzVZcqtYI')

song2 <- list('title' = 'Highway to Hell', 
              'composer' = 'Angus Young, Malcolm Young and Bon Scott', 
              'genre' = 'Hard Rock' , 'year' = 1979, 
              'artist' = 'AC/DC',
              'link' = 'https://www.youtube.com/watch?v=l482T0yNkeo')


# Check the class of song1
print(class(song1))

# We define the class 'song' implicitly. 
attr(song1, 'class') <- 'song'
class(song2) <- 'song'

print(class(song1))

# We now want to define some generic functions for the class 'song'.

# Try methods(print). You will see all classes, for which print() is 
# defined.
methods(print)

# We now define print() for the class 'song'
print.song <- function(mysong) {
    paste(mysong$title, 'is a',mysong$genre,'song by', paste(mysong$artist,'.',sep='')
          , sep = ' ')
}

# This is called a method.

# We can define a new method using UseMethod

get_artist <- function(object) {
  UseMethod('get_artist')
}
get_artist.song <- function(object) {
  return(object$artist)
}

#############################
# 4.2 Exercise              #
#############################

# 1) Write an is.song method.
# 2) Implement a method play. Use the build-in function
#    browseURL.
# 2) Implement a constructor function.


# Possible solution:

# 1)

is.song <- function(mysong) {
  if (class(mysong) == 'song') {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# 2)

play <- function(object) {
  UseMethod('play')
}

play.song <- function(mysong) {
  browseURL(mysong$link)
}

# 3) 

song <- function(title, composer, genre, year, artist) {
  UseMethod('song')
}

song.default <- function(title = NA_character_, composer = NA_character_, 
                 genre = NA_character_ , year = NA_integer_, 
                 artist = NA_character_) {
  
      mysong <- list('title' = title, 'composer' = composer, 
                   'genre' = genre , 'year' = year, 
                   'artist' = artist)
      class(mysong) <- 'song'

      return(mysong)
}

song.pop <- function(title, composer, 
                         genre = 'pop' , year = 2018, 
                         artist) {
  
  mysong <- song.default('genre' = genre , 'year' = year)
  
  return(mysong)
}

  
#############################
# 4.3 formal S4             #
#############################

# Define the class  
setClass(Class = "MyClass" , 
         representation(attr1 = 'numeric',
                        attr2 = 'character'),
         prototype(attr1 = NA_real_, attr2 = NA_character_) 
         )

testobject <- new('MyClass')

# Check out

print(test@attr1)

print(getSlots('MyClass'))

# Define a method

setMethod('methodname',signature = 'classname' , 
          definition = 'function')

#############################
# 4.4 Exercise              #
#############################

# Implement 'song' and play as S4 class

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

play(song1)