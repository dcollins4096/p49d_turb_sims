def write_alignment_file(self, fields, filename):
    '''
    This method writes out fields in the format radmc3d needs to compute
    thermal dust emission. In particular, if you have a field called
    "DustDensity", you can write out a dust_density.inp file.

    Parameters
    ----------

    field : string
    The name of the field to be written out
    filename : string
    The name of the file to write the data to. The filenames radmc3d
    expects for its various modes of operations are described in the
    radmc3d manual.

    '''
    fhandle = open(filename, 'w')

    # write header
    fhandle.write('1 \n')
    fhandle.write(str(self.cell_count) + ' \n')
    if len(field) == 1:
        fhandle.write('1 \n')

    # now write fine layers:
    for layer in self.layers:
        lev = layer.level
        if lev == 0:
            LE = self.domain_left_edge
            N = self.domain_dimensions
        else:
            LE = layer.LeftEdge
            N = layer.ActiveDimensions

    self._write_layer_data_to_file(fhandle, fields, lev, LE, N)

    fhandle.close()

