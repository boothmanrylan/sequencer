import ee


def match_sequence(events, sequence, offset):
    num_events = events.bandNames().length()
    seq_length = sequence.length()

    seq_im = ee.Image.constant(sequence)

    indices = ee.List.sequence(0, num_events.subtract(seq_length))
    indices_im = ee.Image.constant(indices)

    def _check(index, prev):
        compare_to = events.slice(index, seq_length.add(index))
        result = seq_im.eq(compare_to).reduce(ee.Reducer.allNonZero())
        return prev.addBands(result.selfMask())

    checks = ee.Image(indices.iterate(_check, ee.Image())).slice(1)
    matches = checks.multiply(indices_im).add(offset)
    return matches.reduce(ee.Reducer.min())


def match_sequences(events, sequences, offsets):
    indices = ee.List.sequence(0, sequences.length().subtract(1))

    def _match(index):
        curr_seq = sequences.get(index)
        curr_offset = offsets.get(index)
        return match_sequence(events, curr_seq, curr_offset).int()

    matches = ee.ImageCollection.fromImages(indices.map(_match))
    return matches.reduce(ee.Reducer.min())
