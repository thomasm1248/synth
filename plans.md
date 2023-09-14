

# Audio Synth Library

It should have

- A variety of filters and effects
- A variety of wave functions
- A mixer
- An amplitude envelope
- Support for midi?
- Support for plaintext music notation
- Support for full song organization

## Plaintext Music Notation

Ideas

- Partial support for ABC notation

## Instrument development

### Use cases

- Define an instrument that restarts every note being played
- Define an instrument that feeds the previous into an effect chain
- Define an instrument that can lerp between adjacent notes

Words

- Instrument
- Voice
- Note Player
- Track
- Tune
- Note
- Sound

### Voice object

- Constructor sets up oscillators, envelopes, and filters
	- Recieves a note (int where 0 is middle C, +-1 is up/down half-step)
	- Recieves a note duration (float in ms)
- Method for getting the next sample
- Method that returns bool for if it's finished yet
	- Probably delay + note duration + release < samples so far / sample rate * 1000 ms
	- Recommended to not add reverb to voices so that it's easier to know when they're done
- Optional constructor that recieves portamento series instead of single note

### Instrument object

- Has a method that is called each time a new sample needs to be made
	- Uses currently active voices to generate sounds, then passes those through filter objects
- Has a method that is called each time a new note voice is requested by the track schedule
	- Recieves a note (int where 0 is middle C, +-1 is up/down half-step)
	- Recieves a note duration (float in ms)
	- Creates one or more voices

## Synthesis pipeline

### Use cases

- I want to read a wav, apply effects, then write a new wav

### Steps to follow when using the library

- Use library methods to convert strings into instructions for the instruments
- Schedule tunes on the track by calling a library method
	- the measure the tune will start
	- the instrument object to use
	- the tune object (parsed from a string earlier)
- Call the generate method

## Pending Decisions
