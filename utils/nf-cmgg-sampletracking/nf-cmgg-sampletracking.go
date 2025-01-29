package main

import (
	"context"
	"encoding/csv"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"strings"

	smaple_api_client "gitlab.cmgg.be/cmgg/smapleclientgo"

	"github.com/antihax/optional"
	"github.com/google/uuid"
	"github.com/joho/godotenv"
	log "github.com/sirupsen/logrus"
	cli "github.com/urfave/cli/v2"
)

// Structs
type samplesheetSample struct {
	// Sample name
	sample         string `json:"sample"`
	pool           string `json:"pool"`
	sampleBam      string `json:"sample_bam"`
	sampleBamIndex string `json:"sample_bam_index"`
	snpFastq1      string `json:"snp_fastq_1"`
	snpFastq2      string `json:"snp_fastq_2"`
	snpBam         string `json:"snp_bam"`
	snpBamIndex    string `json:"snp_bam_index"`
}

func (s samplesheetSample) header() []string {
	return []string{"sample", "pool", "sample_bam", "sample_bam_index", "snp_fastq_1", "snp_fastq_2", "snp_bam", "snp_bam_index"}
}
func (s samplesheetSample) csv() []string {
	return []string{s.sample, s.pool, s.sampleBam, s.sampleBamIndex, s.snpFastq1, s.snpFastq2, s.snpBam, s.snpBamIndex}
}

func main() {
	// Load environment variables from .env file
	err := godotenv.Load()
	if err != nil {
		log.Warn("Error loading .env file")
	}

	var runName string
	var smapleUsername string
	var smaplePassword string
	var smapleUrl string

	var samplesDir cli.StringSlice
	var snpDir cli.StringSlice

	app := &cli.App{
		Name:  "nf-cmgg-sampletracking",
		Usage: "Generate a sample sheet for the nf-cmgg/sampletracking Nextflow pipeline",
		Flags: []cli.Flag{
			&cli.StringFlag{
				Name:        "run",
				Usage:       "Run Name",
				Aliases:     []string{"r"},
				Required:    true,
				Destination: &runName,
			},
			&cli.StringSliceFlag{
				Name:        "samples_directory",
				Usage:       "Directory containing the full size samples",
				Aliases:     []string{"w"},
				Required:    true,
				Destination: &samplesDir,
			},
			&cli.StringSliceFlag{
				Name:        "snp_directory",
				Usage:       "Directory containing the snp samples",
				Aliases:     []string{"s"},
				Required:    true,
				Destination: &snpDir,
			},
			&cli.StringFlag{
				Name:        "smaple_username",
				Usage:       "Smaple API username",
				EnvVars:     []string{"SMAPLE_USERNAME"},
				Required:    true,
				Destination: &smapleUsername,
			},
			&cli.StringFlag{
				Name:        "smaple_password",
				Usage:       "Smaple API password",
				EnvVars:     []string{"SMAPLE_PASSWORD"},
				Required:    true,
				Destination: &smaplePassword,
			},
			&cli.StringFlag{
				Name:        "smaple_url",
				Usage:       "Smaple Base URL",
				Value:       "https://smaple.cmgg.be",
				EnvVars:     []string{"SMAPLE_URL"},
				Destination: &smapleUrl,
			},
		},
		Action: func(c *cli.Context) error {
			// Authenticate to Smaple
			ctx, client, err := smapleAuthenticate(
				smapleUsername,
				smaplePassword,
				smapleUrl,
			)
			if err != nil {
				log.Fatal("Unable to authenticate to Smaple: ", err)
			}
			// Fetch the samplesheet for <runName>
			samplesheet, resp, err := client.RunApi.GetSampleSheet(ctx, runName)
			if err != nil {
				body, err := io.ReadAll(resp.Body)
				log.Debug("exit code: ", resp.StatusCode, string(body))
				log.Fatal("Failed to get samplesheet for run ", runName, ": ", err)
			}

			// Parse the samplesheet
			var samples = make(map[string]*samplesheetSample)
			for _, lane := range samplesheet.RunLaneObjects {
				for _, pool := range lane.PoolObjects {
					for _, sample := range pool.LibPrepSamplePoolObjects {
						if sample.Tag != "WES" {
							continue
						}
						// Sample_ID
						var sample_id string = ""
						if sample.SampleNumberReference != "" {
							sample_id = sample.SampleNumberReference
						} else {
							sample_id = sample.SampleNumber
						}
						samples[sample_id+"_"+pool.PoolName] = &samplesheetSample{
							sample: sample_id,
							pool:   pool.PoolName,
						}
					}
				}
			}

			// Find all full size (aligned) samples in <samplesDir>
			sampleFiles := []string{}
			for _, dir := range samplesDir.Value() {
				files, err := findFilesByExtension(dir, []string{".bam", ".cram"})
				if err != nil {
					log.Fatal("Error finding files in ", dir, ": ", err)
				}
				sampleFiles = append(sampleFiles, files...)
			}

			// Find all snp samples in <snpDir>
			snpFastq := []string{}
			snpBams := []string{}
			for _, dir := range snpDir.Value() {
				snpfiles, err := findFilesByExtension(dir, []string{".fastq.gz", ".fq.gz"})
				if err != nil {
					log.Fatal("Error finding files in ", dir, ": ", err)
				}
				snpFastq = append(snpFastq, snpfiles...)
				snpbamfiles, err := findFilesByExtension(dir, []string{".bam", ".cram"})
				if err != nil {
					log.Fatal("Error finding files in ", dir, ": ", err)
				}
				snpBams = append(snpBams, snpbamfiles...)
			}

			if err != nil {
				log.Fatal("Error finding files in ", snpDir, ": ", err)
			}

			// Associate each bam/cram with a sample
			for i, sample := range samples {
				// Find the sample bam
				for _, file := range sampleFiles {
					if strings.Contains(file, sample.sample) {
						samples[i].sampleBam = file
						if _, err := os.Stat(file + ".bai"); err == nil {
							samples[i].sampleBamIndex = file + ".bai"
						} else if _, err := os.Stat(file + ".csi"); err == nil {
							samples[i].sampleBamIndex = file + ".csi"
						} else if _, err := os.Stat(file + ".crai"); err == nil {
							samples[i].sampleBamIndex = file + ".crai"
						} else {
							log.Error("Index file not found for ", file)
						}
					}
				}
				// Find the snp fastq
				for _, file := range snpFastq {
					if strings.Contains(file, sample.sample) {
						if strings.Contains(file, "_R1") {
							samples[i].snpFastq1 = file
						} else if strings.Contains(file, "_R2") {
							samples[i].snpFastq2 = file
						}
					}
				}
				// Find the snp bam
				for _, file := range snpBams {
					if strings.Contains(file, sample.sample) {
						samples[i].snpBam = file
						if _, err := os.Stat(file + ".bai"); err == nil {
							samples[i].snpBamIndex = file + ".bai"
						} else if _, err := os.Stat(file + ".csi"); err == nil {
							samples[i].snpBamIndex = file + ".csi"
						} else if _, err := os.Stat(file + ".crai"); err == nil {
							samples[i].snpBamIndex = file + ".crai"
						} else {
							log.Error("Index file not found for ", file)
						}
					}
				}
			}

			// Open the output file
			samplesheet_file, err := os.Create("sampletracking_samplesheet.csv")
			if err != nil {
				log.Fatalf("failed creating file: %s", err)
			}
			samplesheet_writer := csv.NewWriter(samplesheet_file)
			defer samplesheet_writer.Flush()

			samplesheet_writer.Write(samplesheetSample{}.header())
			for _, sample := range samples {
				if sample.sampleBam != "" && ((sample.snpFastq1 != "" && sample.snpFastq2 != "") || sample.snpBam != "") {
					samplesheet_writer.Write(sample.csv())
				}
			}

			return nil
		},
	}

	if err := app.Run(os.Args); err != nil {
		log.Fatal(err)
	}
}

// Authenticate to SMAPLE and return an authenticated client
func smapleAuthenticate(username string, password string, url string) (context.Context, smaple_api_client.APIClient, error) {
	log.Debug("Authenticating to SMAPLE")
	log.Debug("Username:", username)
	log.Debug("URL:", url)

	configuration := smaple_api_client.NewConfiguration()
	configuration.BasePath = url
	client := smaple_api_client.NewAPIClient(configuration)

	opts := &smaple_api_client.AuthApiApiV6AuthLoginPostOpts{
		Body: optional.NewInterface(map[string]string{
			"name":     uuid.New().String(),
			"password": password,
			"user":     username,
		}),
	}

	token, _, err := client.AuthApi.ApiV6AuthLoginPost(context.Background(), opts)
	if err != nil {
		log.Fatal("Unable to fetch Smaple authentication token: ", err)
	}
	log.Debug("Smaple token:", token.AccessToken)

	auth := context.WithValue(context.Background(), smaple_api_client.ContextAPIKey, smaple_api_client.APIKey{
		Key:    token.AccessToken,
		Prefix: "Bearer", // Omit if not necessary.
	})

	return auth, *client, err
}

// Find all files with a given extension in a directory
func findFilesByExtension(dir string, extension []string) ([]string, error) {
	var files []string

	err := filepath.WalkDir(dir, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if !d.IsDir() {
			for _, ext := range extension {
				if strings.HasSuffix(d.Name(), ext) {
					files = append(files, path)
				}
			}
		}
		return nil
	})

	if err != nil {
		return nil, err
	}

	return files, nil
}
