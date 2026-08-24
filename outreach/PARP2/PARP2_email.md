Subject: PARP2 base-editing resistance screen → 3D readout (honest: under-powered, but a clean paralog point)

Hi [name],

I ran the PARP2 arm of the Coelho & Dincer 2024 base-editing resistance screen (olaparib/niraparib) through
BE3D, a structure-function tool from the Iqbal Lab (Broad) and Liau Lab (Harvard) that maps screen
log-fold-change onto 3D structure (AlphaFold Q9UGN5). Independent analysis, and I want to be upfront about
what it can and can't say.

Candidly, the PARP2 run is under-powered. Unlike the companion PARP1 run — where BE3D cleanly re-found the
catalytic pocket zero-shot — PARP2 doesn't discriminate at the residue level: the NAD+/inhibitor pocket isn't
significantly enriched in either direction (only 3 of 9 pocket residues reach significance), precision@10 is
zero, and a no-screen AlphaMissense baseline actually beats BE3D on this protein. The reason is almost
certainly coverage: the PARP2 PARPi arm is ABE-only, 809 guides, so only ~200 residues get scored across a
583-aa protein and the catalytic domain is thinly tiled.

The one solid, useful point is cross-paralog: the PARP2 catalytic pocket is the structural twin of PARP1's
(H428/Y462/E558 ↔ H862/Y896/E988), and in both paralogs it sits in the depletion/essential direction, never
the resistance direction — mutating the pocket is deleterious, not protective, in both. That conserved
signature is real; BE3D just recovers it strongly in PARP1 and weakly in PARP2.

Two practical notes if you revisit this: (1) the screen numbers PARP2 in a −13 short-isoform frame, so map
carefully to canonical/structure; (2) a CBE arm plus more guides would likely make PARP2 as informative as
PARP1.

Attached: the flagged hotspot table (with the under-power caveat), G2P files for 3D viewing, and a one-page
brief. Happy to re-run if a deeper PARP2 screen exists.

Best,
[name]
