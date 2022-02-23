"""Methods that find earliest occurences of exact sequences in images."""

import ee


def match_sequence(events, sequence, offset=0):
    """Finds the earliest occurence of sequence in events.

    Determines which, if any, pixels contain sequence in them and at what
    band index in events the sequence first appears (plus offset to allow
    the point of interest to occur part way through the sequence).

    Args:
        events: An ee.Image whose bands represent classifications through time.
        sequence: An ee.List of class values to search for e.g. [1, 2, 2, 2].
        offset: An integer representing the offset of from the start of the
            sequence to the point of interest in the sequence. Default to 0 to
            return the start of the sequence.

    Returns:
        An ee.Image with one band whose values are the band index of the input
        image where sequence first occurs plus offset
    """
    events = ee.Image(events)
    sequence = ee.List(sequence)
    offset = ee.Number(offset)

    num_events = events.bandNames().length()
    seq_length = sequence.length()

    seq_im = ee.Image.constant(sequence)

    indices = ee.List.sequence(0, num_events.subtract(seq_length))
    indices_im = ee.Image.constant(indices)

    def _check(index, prev):
        prev = ee.Image(prev)
        compare_to = events.slice(index, seq_length.add(index))
        result = seq_im.eq(compare_to).reduce(ee.Reducer.allNonZero())
        return prev.addBands(result.selfMask())

    # slice(1) to drop the blank image iterate starts with
    checks = ee.Image(indices.iterate(_check, ee.Image())).slice(1)
    matches = checks.multiply(indices_im).add(offset)
    return matches.reduce(ee.Reducer.min())


def match_sequences(events, sequences, offsets):
    """Finds the earliest occurence of any of sequences in events.

    Determines which, if any, pixels in events contain one of sequences in
    them and at what band index in events the first of sequences appears (plus
    some offset to allow the point of interest to occur part way through a
    sequence).

    Args:
        events: An ee.Image whose bands represent classifications through time.
        sequences: An ee.List of ee.List of integers representing each of the
            sequences of class values to search for.
        offsets: An ee.List of intergers representing the offsets for each of
            sequences.

    Returns:
        An ee.Image with one band whose values represent the band index of
        events where the first of sequences occurs plus some sequence specific
        offset.
    """
    events = ee.Image(events)
    sequences = ee.List(sequences)
    offsets = ee.List(offsets)

    indices = ee.List.sequence(0, sequences.length().subtract(1))

    def _match(index):
        curr_seq = sequences.get(index)
        curr_offset = offsets.get(index)
        return match_sequence(events, curr_seq, curr_offset).int()

    matches = ee.ImageCollection.fromImages(indices.map(_match))
    return matches.reduce(ee.Reducer.min())
