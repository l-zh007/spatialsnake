














      sq.gr.co_occurrence(adata, cluster_key="clusters")
      sq.pl.co_occurrence(
        adata,
        cluster_key="clusters",
        clusters="0",
        figsize=(8, 4),
        )

      plt.savefig(
          os.path.join(output_dir, f"{key}co_occur.png"),
          dpi=300,
          bbox_inches='tight')
      plt.show()
      plt.close()
      
      
      sq.gr.ligrec(
      adata,
      n_perms=100,
      cluster_key="clusters",
      use_raw=False)
      
      
      
      
      sq.pl.ligrec(
        adata,
        cluster_key="clusters",
        source_groups="0",
        target_groups=["1", "2"],
        means_range=(3, np.inf),
        alpha=1e-4,
        swap_axes=True)


      plt.savefig(
          os.path.join(output_dir, f"{key}co_occur2.png"),
          dpi=300,
          bbox_inches='tight')
      plt.show()
      plt.close()
      
      
      
      
      
      
      
      

    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key="clusters")
    sq.pl.nhood_enrichment(adata, cluster_key="clusters", figsize=(5, 5))
    plt.savefig(
        os.path.join(output_dir, f"{key}Cells_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()

    mode = "L"
    sq.gr.ripley(adata, cluster_key="clusters", mode=mode, max_dist=500)
    sq.pl.ripley(adata, cluster_key="clusters", mode=mode)
    plt.savefig(
        os.path.join(output_dir, f"{key}Across_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()
