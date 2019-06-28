#
# Show the hidden part and hide the shown part of a surface.
# Can be used with the hide dust tool to show only the dust to mask density
# maps to set the grid points in the dust blobs to zero.
#
import Surface
plist = Surface.selected_surface_pieces()
for p in plist:
    m = p.triangleAndEdgeMask
    if m is not None:
        minv = m^8
        p.triangleAndEdgeMask = minv
        Surface.set_visibility_method('invert', p.model) # Suppress surface mask auto updates
