#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>
#include <vector>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "caffe/data_transformer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/layers/ndim_data_layer.hpp"
#include "caffe/util/benchmark.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/util/rng.hpp"



namespace caffe {

template <typename Dtype>
NDimDataLayer<Dtype>::~NDimDataLayer<Dtype>() {
  this->StopInternalThread();
}

//create simple shape vector from blobshape
template <typename Dtype>
const vector<int>  NDimDataLayer<Dtype>::blob2vec(const BlobShape& b) const
{
  CHECK_LE(b.dim_size(), kMaxBlobAxes);
  vector<int> shape_vec(b.dim_size());
  for (int i = 0, n = b.dim_size(); i < n; ++i) {
    shape_vec[i] = b.dim(i);
  }
  return shape_vec;
}

template <typename Dtype>
void NDimDataLayer<Dtype>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  string root_folder = this->layer_param_.ndim_data_param().root_folder();
  const bool balanced  = this->layer_param_.ndim_data_param().balanced();
  num_rotations = this->layer_param_.ndim_data_param().rotate();
  all_pos_ = actives_pos_ = decoys_pos_ = 0;

  //shape must come from parameters
  const int batch_size = this->layer_param_.ndim_data_param().batch_size();
  CHECK_GT(batch_size, 0) << "Positive batch size required";

  if(!this->layer_param_.ndim_data_param().inmemory())
  {
		// Read the file with filenames and labels
		const string& source = this->layer_param_.ndim_data_param().source();

		LOG(INFO) << "Opening file " << source;
		std::ifstream infile(source.c_str());
		CHECK((bool)infile) << "Could not load " << source;


		string line, fname;
		vector<std::string> binmaps;

		while (getline(infile, line)) {
			stringstream example(line);
			float label = 0;
			//first the label
			example >> label;
			//then all binmaps for the example
			binmaps.clear();

			while(example >> fname) {
				if(fname.length() > 0 && fname[0] == '#')
					break; //ignore rest of line
				binmaps.push_back(fname);
			}

			if(binmaps.size() == 0) //ignore empty lines
				continue;

			all_.push_back(make_pair(binmaps,label));
			if(label) actives_.push_back(binmaps);
			else decoys_.push_back(binmaps);

			if(label != 0.0 && label != 1.0) {
			  CHECK(!balanced) << "Non-binary labels with balanced set to true";
			}
		}

		if (this->layer_param_.ndim_data_param().shuffle()) {
			// randomly shuffle data
			LOG(INFO) << "Shuffling data";
			const unsigned int prefetch_rng_seed = caffe_rng_rand();
			prefetch_rng_.reset(new Caffe::RNG(prefetch_rng_seed));
			Shuffle();
		}
		LOG(INFO) << "A total of " << all_.size() << " examples.";

		// Check if we would need to randomly skip a few data points
		if (this->layer_param_.ndim_data_param().rand_skip()) {
			unsigned int skip = caffe_rng_rand() %
					this->layer_param_.ndim_data_param().rand_skip();
			LOG(INFO) << "Skipping first " << skip << " data points.";
			CHECK_GT(all_.size(), skip) << "Not enough points to skip";
			all_pos_ = skip;
			actives_pos_ = skip % actives_.size();
			decoys_pos_ = skip % decoys_.size();
		}

		if(balanced) {
			CHECK_GT(batch_size, 1) << "Batch size must be > 1 with balanced option.";
		}
  }

  vector<int> example_shape = blob2vec(this->layer_param_.ndim_data_param().shape());
  top_shape.clear();
  top_shape.push_back(1);

  example_size =1;
  for(unsigned i = 0, n = example_shape.size(); i < n; i++) {
    CHECK_GT(example_shape[i], 0) << "Positive shape dimension required";
    top_shape.push_back(example_shape[i]);
    example_size *= example_shape[i];
  }
  //shape of single data
  this->transformed_data_.Reshape(top_shape);

  // Reshape prefetch_data and top[0] according to the batch_size.
  top_shape[0] = batch_size;
  for (int i = 0; i < this->prefetch_.size(); ++i) {
    this->prefetch_[i]->data_.Reshape(top_shape);
  }
  top[0]->Reshape(top_shape);

  // label
  vector<int> label_shape(1, batch_size);
  top[1]->Reshape(label_shape);
  for (int i = 0; i < this->prefetch_.size(); ++i) {
    this->prefetch_[i]->label_.Reshape(label_shape);
  }
}

template <typename Dtype>
void NDimDataLayer<Dtype>::Shuffle() {
  caffe::rng_t* prefetch_rng =
      static_cast<caffe::rng_t*>(prefetch_rng_->generator());
    shuffle(actives_.begin(), actives_.end(), prefetch_rng);
    shuffle(decoys_.begin(), decoys_.end(), prefetch_rng);
    shuffle(all_.begin(), all_.end(), prefetch_rng);
}

//copy raw floating point data from files into buffer
//files should represent a single example
template <typename Dtype>
void  NDimDataLayer<Dtype>::load_data_from_files(Dtype* buffer, const std::string& root, const vector<std::string>& files)
{
  using namespace boost::iostreams;

  CHECK_GT(files.size(), 0) << "Missing binmaps files";
  
  unsigned long total = 0;
  for(unsigned i = 0, n = files.size(); i < n; i++)
  {
    std::string fname = root + files[i];

    std::ifstream inorig(fname.c_str(), std::ios::binary);
    CHECK((bool)inorig) << "Could not load " << fname;
    std::string::size_type pos = fname.rfind(".gz");
    if (pos != std::string::npos)
    {
      //setup streams to do copy
      basic_array_sink<char> data(((char*)buffer)+total, example_size*sizeof(Dtype)-total);
      filtering_stream<input> in;
      in.push(gzip_decompressor());
      in.push(inorig);
      total += copy(in, data);

      CHECK_EQ(inorig.peek(), EOF) << "File " << fname << " not fully read.  Are grid input sizes correct?";
    }
    else
    {
      //do direct read of full file
      //get length of file
      inorig.seekg(0, inorig.end);
      std::streamsize size = inorig.tellg();
      inorig.seekg(0, inorig.beg);
      CHECK_LE(total+size, example_size*sizeof(Dtype)) << "file " << fname << " is too big";
      char *data = ((char*)buffer)+total; //starting place to copy data
      inorig.read(data, size);
      total += size;
    }

  }

  CHECK_EQ(total,example_size*sizeof(Dtype)) << "Incorrect size of inputs (" << total << " vs. " << example_size*sizeof(Dtype) << ") on " << files[0];

  if(current_rotation > 0) {
    rotate_data(buffer, current_rotation);
  }
  if( this->layer_param_.ndim_data_param().check()) {
    for(unsigned i = 0; i < example_size; i++) {
      CHECK(finite(buffer[i])) << "Not finite value at " << i << " in " << files[0];
    }

  }
}

static unsigned rotate_coords(unsigned i, unsigned j, unsigned k, const unsigned I, const unsigned J, const unsigned K,unsigned rot)
{
  CHECK_LT(rot,24) << "Invalid rotation " << rot << " (must be <24)";
  CHECK_EQ(I,J) << "Require cube input for rotate"; //could remove this requirement by adjust IJK appropriately
  CHECK_EQ(J,K) << "Require cube input for rotate";
  unsigned newi = i, newj = j, newk = k;

  //rotate to a face
  switch(rot%6) {
  case 0:
    newi = i; newj = j; newk = k;
    break;
  case 1:
    newi = j;
    newj = I-i-1;
    newk = k;
    break;
  case 2:
    newi = I-i-1;
    newj = J-j-1;
    newk = k;
    break;
  case 3:
    newi = J-j-1;
    newj = i;
    newk = k;
    break;
  case 4:
    newi = k;
    newk = I-i-1;
    newj = j;
    break;
  case 5:
    newi = K-k-1;
    newk = i;
    newj = j;
    break;
  }
  i = newi;
  j = newj;
  k = newk;

  //now rotate around other axis
  rot /= 6;

  switch(rot%4) {
  case 0:
    newi = i; newj = j; newk = k;
    break;
  case 1:
    newj = k;
    newk = J-j-1;
    newi = i;
    break;
  case 2:
    newj = J-j-1;
    newk = K-k-1;
    newi = i;
    break;
  case 3:
    newj = K-k-1;
    newk = j;
    newi = i;
    break;
  }
  i = newi;
  j = newj;
  k = newk;

  return ((newi*J)+newj)*K+k;
}

// This takes a buffer of a single example (example_sie elements)
//and will perform a rotation along the axis; rot specifies which of the 32
//possible rotations will be performed
template <typename Dtype>
void NDimDataLayer<Dtype>::rotate_data(Dtype *data, unsigned rot) {
  Dtype *tdata = (Dtype*)malloc(sizeof(Dtype)*example_size); //transformed data

  CHECK_EQ(top_shape.size(),5) << "Rotate only available for 3D inputs";
  const unsigned I = top_shape[2];
  const unsigned J = top_shape[3];
  const unsigned K = top_shape[4];

  Dtype *from = data;
  Dtype *to = tdata;

  for(unsigned c = 0, nc = top_shape[1]; c < nc; c++) {
    //for each channel
    for(unsigned i = 0; i < I; i++) {
      for(unsigned j = 0; j < J; j++) {
        for(unsigned k = 0; k < K; k++) {
          unsigned origpos =  ((i*J)+j)*K+k;
          unsigned newpos = rotate_coords(i,j,k,I,J,K,rot);
          CHECK_LT(newpos+c*I*J*K,example_size) << "out of bounds " << c <<","<<i<<"," <<j<<","<<k<<" rot " << rot <<"\n";
          to[newpos] = from[origpos];
        }
      }
    }

    //next channel
    from += I*J*K;
    to += I*J*K;
  }

  memcpy(data, tdata,example_size*sizeof(Dtype));
  free(tdata);
}


template <typename Dtype>
void NDimDataLayer<Dtype>::memoryIsSet()
{
	boost::unique_lock<boost::mutex> lock(mem_mutex);
	unsigned batch_size = top_shape[0];
	unsigned add = 1;

	if(num_rotations > 0) {
		CHECK_LE(batch_size, num_rotations);
		CHECK_EQ(num_rotations % batch_size, 0);
		add = num_rotations/batch_size;
	}
	data_avail += add;

	mem_cond.notify_one();
}

// This function is called on prefetch thread
template <typename Dtype>
void NDimDataLayer<Dtype>::load_batch(Batch<Dtype>* batch) {
  CPUTimer batch_timer;
  batch_timer.Start();
  string root_folder = this->layer_param_.ndim_data_param().root_folder();
  const bool balanced  = this->layer_param_.ndim_data_param().balanced();
  const bool inmem = this->layer_param_.ndim_data_param().inmemory();

  CHECK(batch->data_.count());
  CHECK(this->transformed_data_.count());
  unsigned batch_size = top_shape[0];
  CHECK_GT(batch_size, 0) << "Positive batch size required";
  vector<int> offind(1, 0);

  batch->data_.Reshape(top_shape);

  Dtype* prefetch_data = batch->data_.mutable_cpu_data();
  Dtype* prefetch_label = batch->label_.mutable_cpu_data();

  if(inmem) {
  	boost::unique_lock<boost::mutex> lock(mem_mutex);
  	while(data_avail == 0)
  	{
  		mem_cond.wait(lock);
  	}
  	data_avail--;

  	//memory is now available
  	CHECK_EQ(batch->data_.count() % memdata.size() , 0) << "Mismatch in ndim memory size";
  	unsigned batchsz = batch->data_.count() / memdata.size();

  	for(unsigned b = 0; b < batchsz; b++) {
  		Dtype *outdata = prefetch_data+b*memdata.size();
			memcpy(outdata, &memdata[0], memdata.size()*sizeof(float));

			if(current_rotation > 0) {
				rotate_data(outdata, current_rotation);
			}

			if (num_rotations > 0) {
				current_rotation = (current_rotation+1)%num_rotations;
			}
  	}
  }
  else {
		if(balanced) { //load equally from actives/decoys
			unsigned nactives = batch_size/2;

			int item_id = 0;
			unsigned asz = actives_.size();
			for (item_id = 0; item_id < nactives; ++item_id) {
				offind[0] = item_id;
				int offset = batch->data_.offset(offind);
				load_data_from_files(prefetch_data+offset, root_folder, actives_[actives_pos_]);
				prefetch_label[item_id] = 1;

				actives_pos_++;
				if(actives_pos_ >= asz) {
					DLOG(INFO) << "Restarting actives data prefetching from start.";
					actives_pos_ = 0;
					if (this->layer_param_.ndim_data_param().shuffle()) {
						shuffle(actives_.begin(), actives_.end(), static_cast<caffe::rng_t*>(prefetch_rng_->generator()));
					}
					if (num_rotations > 0) {
						current_rotation = (current_rotation+1)%num_rotations;
					}
				}
			}
			unsigned dsz = decoys_.size();
			for (; item_id < batch_size; ++item_id) {
				offind[0] = item_id;
				int offset = batch->data_.offset(offind);
				load_data_from_files(prefetch_data+offset, root_folder, decoys_[decoys_pos_]);
				prefetch_label[item_id] = 0;

				decoys_pos_++;
				if(decoys_pos_ >= dsz) {
					DLOG(INFO) << "Restarting decoys data prefetching from start.";
					decoys_pos_ = 0;
					if (this->layer_param_.ndim_data_param().shuffle()) {
						shuffle(decoys_.begin(), decoys_.end(), static_cast<caffe::rng_t*>(prefetch_rng_->generator()));
					}
					if (num_rotations > 0) {
						current_rotation = (current_rotation+1)%num_rotations;
					}
				}
			}

		} else {
			//load from all
			unsigned sz = all_.size();
			for (int item_id = 0; item_id < batch_size; ++item_id) {
				offind[0] = item_id;
				int offset = batch->data_.offset(offind);
				load_data_from_files(prefetch_data+offset, root_folder, all_[all_pos_].first);
				prefetch_label[item_id] = all_[all_pos_].second;

				all_pos_++;
				if(all_pos_ >= sz) {
					DLOG(INFO) << "Restarting data prefetching from start.";
					all_pos_ = 0;
					if (this->layer_param_.ndim_data_param().shuffle()) {
						shuffle(all_.begin(), all_.end(), static_cast<caffe::rng_t*>(prefetch_rng_->generator()));
					}
					if (num_rotations > 0) {
						current_rotation = (current_rotation+1)%num_rotations;
					}
				}
			}
		}
  }

  batch_timer.Stop();
  DLOG(INFO) << "Prefetch batch: " << batch_timer.MilliSeconds() << " ms.";
}

INSTANTIATE_CLASS(NDimDataLayer);
REGISTER_LAYER_CLASS(NDimData);

}  // namespace caffe
